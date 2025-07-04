/// Liftover support for alt-aware mapping
/// 
/// This module implements DRAGMAP-compatible liftover functionality
/// for mapping between primary and alternate reference contigs.

use std::collections::HashMap;
use std::path::Path;
use anyhow::{Result, anyhow};
use serde::{Deserialize, Serialize};

/// Liftover status codes (matches DRAGMAP's implementation)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum LiftoverStatus {
    /// Position is unaligned
    Unaligned = 0,
    /// Position matches between alt and primary
    Match = 1,
    /// Position is an insertion in alt relative to primary
    Insert = 2,
    /// Position has mixed status (complex structural variation)
    Mixed = 3,
}

/// Liftover direction for coordinate transformation
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LiftoverDirection {
    /// Forward direction (positive strand)
    Forward = 1,
    /// Reverse direction (negative strand)
    Reverse = -1,
}

/// Liftover node in the hierarchical trie structure
/// Based on DRAGMAP's liftoverNode_t structure
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LiftoverNode {
    /// Position in the reference
    pub position: u64,
    /// Width of the region covered by this node
    pub width: u16,
    /// Direction (+1 for forward, -1 for reverse)
    pub direction: i8,
    /// Status of this liftover region
    pub status: LiftoverStatus,
    /// Child nodes for hierarchical structure
    pub children: Vec<LiftoverNode>,
}

/// Liftover group representing primary/alt contig relationships
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LiftoverGroup {
    /// Unique group identifier (28-bit value)
    pub group_id: u32,
    /// Primary contig information
    pub primary_contig: ContigInfo,
    /// Alternative contigs in this group
    pub alt_contigs: Vec<ContigInfo>,
    /// Liftover mappings from alt to primary coordinates
    pub liftover_mappings: HashMap<String, Vec<LiftoverNode>>,
}

/// Information about a reference contig
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContigInfo {
    /// Contig name (e.g., "chr1", "chr1_KI270706v1_random")
    pub name: String,
    /// Length of the contig in bases
    pub length: u64,
    /// Whether this is an alternate contig
    pub is_alt: bool,
    /// Flat reference start position (0-based)
    pub flat_start: u64,
    /// Flat reference end position (0-based, exclusive)
    pub flat_end: u64,
}

/// Liftover code for hash table records (matches DRAGMAP's ExtendTableRecord)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LiftCode {
    /// Regular primary contig
    None = 0,
    /// Alternate contig
    Alt = 1,
    /// Primary contig
    Pri = 2,
    /// Different primary contig
    DifPri = 3,
}

/// Main liftover manager for alt-aware mapping
#[derive(Debug)]
pub struct LiftoverManager {
    /// Map from contig name to liftover group
    groups: HashMap<String, LiftoverGroup>,
    /// Threshold seed index beginning alternate contigs
    alt_seed_threshold: u32,
    /// Threshold flat reference position beginning alternate contigs
    alt_start_threshold: u64,
    /// Map from flat position to contig info
    position_to_contig: Vec<(u64, ContigInfo)>,
    /// Map from group ID to group
    group_by_id: HashMap<u32, LiftoverGroup>,
}

impl LiftoverNode {
    /// Create a new liftover node
    pub fn new(position: u64, width: u16, direction: i8, status: LiftoverStatus) -> Self {
        Self {
            position,
            width,
            direction,
            status,
            children: Vec::new(),
        }
    }
    
    /// Add a child node to this liftover node
    pub fn add_child(&mut self, child: LiftoverNode) {
        self.children.push(child);
    }
    
    /// Find the liftover node covering a specific position
    pub fn find_covering_node(&self, target_position: u64) -> Option<&LiftoverNode> {
        // Check if this node covers the target position
        if target_position >= self.position && target_position < self.position + self.width as u64 {
            // Check children for more specific coverage
            for child in &self.children {
                if let Some(covering_child) = child.find_covering_node(target_position) {
                    return Some(covering_child);
                }
            }
            // No child covers it, so this node is the covering node
            return Some(self);
        }
        None
    }
    
    /// Liftover a position using this node's mapping
    pub fn liftover_position(&self, alt_position: u64) -> Option<(u64, bool)> {
        if let Some(_covering_node) = self.find_covering_node(alt_position) {
            match self.status {
                LiftoverStatus::Match => {
                    // Direct coordinate mapping
                    let offset = alt_position - self.position;
                    let primary_position = if self.direction > 0 {
                        self.position + offset
                    } else {
                        self.position + self.width as u64 - offset - 1
                    };
                    Some((primary_position, self.direction < 0))
                }
                LiftoverStatus::Insert => {
                    // Position is an insertion in alt - map to closest primary position
                    Some((self.position, self.direction < 0))
                }
                LiftoverStatus::Unaligned | LiftoverStatus::Mixed => {
                    // Cannot reliably liftover
                    None
                }
            }
        } else {
            None
        }
    }
}

impl LiftoverGroup {
    /// Create a new liftover group
    pub fn new(group_id: u32, primary_contig: ContigInfo) -> Self {
        Self {
            group_id,
            primary_contig,
            alt_contigs: Vec::new(),
            liftover_mappings: HashMap::new(),
        }
    }
    
    /// Add an alternative contig to this group
    pub fn add_alt_contig(&mut self, alt_contig: ContigInfo, liftover_nodes: Vec<LiftoverNode>) {
        let contig_name = alt_contig.name.clone();
        self.alt_contigs.push(alt_contig);
        self.liftover_mappings.insert(contig_name, liftover_nodes);
    }
    
    /// Check if a contig belongs to this group
    pub fn contains_contig(&self, contig_name: &str) -> bool {
        if self.primary_contig.name == contig_name {
            return true;
        }
        self.alt_contigs.iter().any(|alt| alt.name == contig_name)
    }
    
    /// Liftover position from alt contig to primary coordinates
    pub fn liftover_to_primary(&self, alt_contig: &str, alt_position: u64) -> Option<(u64, bool)> {
        if let Some(liftover_nodes) = self.liftover_mappings.get(alt_contig) {
            // Find the appropriate liftover node for this position
            for node in liftover_nodes {
                if let Some(result) = node.liftover_position(alt_position) {
                    return Some(result);
                }
            }
        }
        None
    }
    
    /// Get all contigs in this group (primary + alts)
    pub fn get_all_contigs(&self) -> Vec<&ContigInfo> {
        let mut contigs = vec![&self.primary_contig];
        contigs.extend(self.alt_contigs.iter());
        contigs
    }
}

impl LiftoverManager {
    /// Create a new liftover manager
    pub fn new() -> Self {
        Self {
            groups: HashMap::new(),
            alt_seed_threshold: 0,
            alt_start_threshold: 0,
            position_to_contig: Vec::new(),
            group_by_id: HashMap::new(),
        }
    }
    
    /// Load liftover configuration from SAM format liftover file
    pub fn load_from_sam<P: AsRef<Path>>(path: P) -> Result<Self> {
        let _path_ref = path.as_ref();
        
        // TODO: Parse actual SAM liftover file
        // For now, create a minimal liftover manager for testing
        
        let mut manager = Self::new();
        
        // Create a simple example liftover group for chr1
        let primary_contig = ContigInfo {
            name: "chr1".to_string(),
            length: 248956422,
            is_alt: false,
            flat_start: 0,
            flat_end: 248956422,
        };
        
        let mut group = LiftoverGroup::new(1, primary_contig);
        
        // Add an example alt contig
        let alt_contig = ContigInfo {
            name: "chr1_KI270706v1_random".to_string(),
            length: 175055,
            is_alt: true,
            flat_start: 248956422,
            flat_end: 248956422 + 175055,
        };
        
        // Create simple liftover nodes for this alt contig (split into chunks due to u16 width limit)
        let liftover_nodes = vec![
            LiftoverNode::new(0, 65535, 1, LiftoverStatus::Match),
            LiftoverNode::new(65535, 65535, 1, LiftoverStatus::Match),
            LiftoverNode::new(131070, 43985, 1, LiftoverStatus::Match),
        ];
        
        group.add_alt_contig(alt_contig, liftover_nodes);
        
        manager.add_group(group)?;
        manager.alt_start_threshold = 248956422;
        
        Ok(manager)
    }
    
    /// Add a liftover group to the manager
    pub fn add_group(&mut self, group: LiftoverGroup) -> Result<()> {
        let group_id = group.group_id;
        
        // Add primary contig
        self.groups.insert(group.primary_contig.name.clone(), group.clone());
        
        // Add primary contig to position mapping
        self.position_to_contig.push((group.primary_contig.flat_start, group.primary_contig.clone()));
        
        // Add alt contigs
        for alt_contig in &group.alt_contigs {
            self.groups.insert(alt_contig.name.clone(), group.clone());
            
            // Add to position mapping
            self.position_to_contig.push((alt_contig.flat_start, alt_contig.clone()));
        }
        
        // Add to group lookup
        self.group_by_id.insert(group_id, group);
        
        // Sort position mapping for binary search
        self.position_to_contig.sort_by_key(|(pos, _)| *pos);
        
        Ok(())
    }
    
    /// Get liftover group for a contig
    pub fn get_group_for_contig(&self, contig_name: &str) -> Option<&LiftoverGroup> {
        self.groups.get(contig_name)
    }
    
    /// Get liftover group by ID
    pub fn get_group_by_id(&self, group_id: u32) -> Option<&LiftoverGroup> {
        self.group_by_id.get(&group_id)
    }
    
    /// Check if a position is in an alternate contig
    pub fn is_alt_position(&self, flat_position: u64) -> bool {
        flat_position >= self.alt_start_threshold
    }
    
    /// Get contig info for a flat reference position
    pub fn get_contig_for_position(&self, flat_position: u64) -> Option<&ContigInfo> {
        // Binary search for the contig containing this position
        let idx = self.position_to_contig
            .binary_search_by(|(pos, _)| pos.cmp(&flat_position))
            .unwrap_or_else(|i| i.saturating_sub(1));
            
        if let Some((_, contig)) = self.position_to_contig.get(idx) {
            if flat_position >= contig.flat_start && flat_position < contig.flat_end {
                return Some(contig);
            }
        }
        
        None
    }
    
    /// Get liftover code for a flat reference position
    pub fn get_lift_code(&self, flat_position: u64) -> LiftCode {
        if let Some(contig) = self.get_contig_for_position(flat_position) {
            if contig.is_alt {
                LiftCode::Alt
            } else {
                LiftCode::Pri
            }
        } else {
            LiftCode::None
        }
    }
    
    /// Liftover position from alt to primary coordinates
    pub fn liftover_position(&self, contig_name: &str, position: u64) -> Option<(String, u64, bool)> {
        if let Some(group) = self.get_group_for_contig(contig_name) {
            if let Some((primary_pos, is_reverse)) = group.liftover_to_primary(contig_name, position) {
                return Some((group.primary_contig.name.clone(), primary_pos, is_reverse));
            }
        }
        None
    }
    
    /// Get all liftover groups
    pub fn get_all_groups(&self) -> Vec<&LiftoverGroup> {
        self.group_by_id.values().collect()
    }
    
    /// Set alt contig thresholds
    pub fn set_alt_thresholds(&mut self, seed_threshold: u32, position_threshold: u64) {
        self.alt_seed_threshold = seed_threshold;
        self.alt_start_threshold = position_threshold;
    }
}

impl Default for LiftoverManager {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_liftover_node() {
        let mut node = LiftoverNode::new(1000, 500, 1, LiftoverStatus::Match);
        assert_eq!(node.position, 1000);
        assert_eq!(node.width, 500);
        assert_eq!(node.direction, 1);
        
        // Test position coverage
        assert!(node.find_covering_node(1200).is_some());
        assert!(node.find_covering_node(900).is_none());
        assert!(node.find_covering_node(1500).is_none());
        
        // Test liftover
        let (lifted_pos, is_reverse) = node.liftover_position(1200).unwrap();
        assert_eq!(lifted_pos, 1200);
        assert!(!is_reverse);
    }
    
    #[test]
    fn test_liftover_group() {
        let primary = ContigInfo {
            name: "chr1".to_string(),
            length: 1000000,
            is_alt: false,
            flat_start: 0,
            flat_end: 1000000,
        };
        
        let mut group = LiftoverGroup::new(1, primary);
        assert_eq!(group.group_id, 1);
        assert_eq!(group.primary_contig.name, "chr1");
        
        let alt = ContigInfo {
            name: "chr1_alt".to_string(),
            length: 50000,
            is_alt: true,
            flat_start: 1000000,
            flat_end: 1050000,
        };
        
        let liftover_nodes = vec![
            LiftoverNode::new(0, 50000, 1, LiftoverStatus::Match),
        ];
        
        group.add_alt_contig(alt, liftover_nodes);
        assert!(group.contains_contig("chr1"));
        assert!(group.contains_contig("chr1_alt"));
        assert!(!group.contains_contig("chr2"));
    }
    
    #[test]
    fn test_liftover_manager() {
        let mut manager = LiftoverManager::new();
        
        let primary = ContigInfo {
            name: "chr1".to_string(),
            length: 1000000,
            is_alt: false,
            flat_start: 0,
            flat_end: 1000000,
        };
        
        let group = LiftoverGroup::new(1, primary);
        manager.add_group(group).unwrap();
        
        assert!(manager.get_group_for_contig("chr1").is_some());
        assert!(manager.get_group_by_id(1).is_some());
        assert_eq!(manager.get_lift_code(500000), LiftCode::Pri);
    }
}