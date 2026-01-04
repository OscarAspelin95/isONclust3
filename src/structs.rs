use std::collections::HashSet;
use std::fmt;

pub enum Cluster<T, U> {
    read_ids(HashSet<T>),
    mini_seqs(HashSet<U>),
}

/// Represents a minimizer along with its starting position in the input string.
/// TODO: rename to indexer or similar
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Minimizer {
    pub sequence: String,
    pub position: usize,
    // pub is_representative: bool
}

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Minimizer_hashed {
    pub sequence: u64,
    pub position: usize,
    // pub is_representative: bool
}

pub(crate) struct GtfEntry {
    pub seqname: String,
    pub source: String,
    pub feature: String,
    pub start: usize,
    pub end: usize,
    //pub score: f64,
    pub strand: bool,
    //pub frame: i8,
    //pub attribute: String
}
impl fmt::Display for GtfEntry {
    // enables displaying the fasta record
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}, {}, {}", self.seqname, self.source, self.feature)
    }
}

#[derive(Debug, Eq, PartialEq, Hash, Clone)]
pub(crate) struct Coord_obj {
    pub startpos: u64,
    pub endpos: u64,
}
impl fmt::Display for Coord_obj {
    // enables displaying the fasta record
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}, {}", self.startpos, self.endpos)
    }
}
impl Coord_obj {
    pub fn overlaps_with(&self, other: &Coord_obj) -> bool {
        let mut overlaps = false;
        if (self.startpos <= other.endpos && self.endpos >= other.startpos)
            || (other.startpos <= self.endpos && other.endpos >= self.startpos)
        {
            overlaps = true;
        }
        overlaps
    }
}

#[derive(Clone)]
pub(crate) struct FastaRecord {
    //a struct used to store fasta records
    pub header: String,
    pub sequence: String,
}

impl fmt::Display for FastaRecord {
    // enables displaying the fasta record
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}\n{}", self.header, self.sequence)
    }
}

pub(crate) struct FastqRecord {
    //a struct used to store fastq records
    header: String,
    sequence: String,
    //quality_header: String,
    //quality: String,
}
impl fmt::Display for FastqRecord {
    // enables displaying the fastq record
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}\n{}", self.header, self.sequence)
    }
}

pub(crate) struct FastqRecord_isoncl_init {
    //a struct used to store fastq records
    pub header: String,
    pub sequence: String,
    pub quality: String,
    pub score: f64,
    pub error_rate: f64,
}
pub(crate) struct internal_gff {
    pub seqname: String,
    pub feature_type: String,
    pub start_coord: u64,
    pub end_coord: u64,
}

impl FastqRecord_isoncl_init {
    pub fn get_header(&self) -> &str {
        &self.header
    }
    //pub fn get_int_id(&self)->&i32{
    // &self.internal_id
    //}
    pub fn get_sequence(&self) -> &str {
        &self.sequence
    }
    pub fn get_quality(&self) -> &str {
        &self.quality
    }
    //pub fn get_score(&self)->&f64{ &self.score }
    pub fn get_err_rate(&self) -> &f64 {
        &self.error_rate
    }
    pub fn set_error_rate(&mut self, new_error_rate: f64) {
        self.error_rate = new_error_rate
    }
    pub fn set_score(&mut self, new_score: f64) {
        self.score = new_score
    }
}

impl fmt::Display for FastqRecord_isoncl_init {
    // enables displaying the fastq record
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}\n{}", self.header, self.sequence)
    }
}
