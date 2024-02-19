use std::collections::VecDeque;
use std::ops::Index;
use rayon::prelude::*;
//use crate::file_actions::FastqRecord_isoncl_init;
use std::cmp::max;
use crate::structs::{FastqRecord_isoncl_init, FastaRecord, Minimizer, Minimizer_hashed};
use crate::clustering::reverse_complement;
use std::borrow::{Borrow, Cow};
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};



//takes an object T and hashes it via DefaultHasher. Used to improve search for minimizers in the data
pub fn calculate_hash<T: Hash>(t: &T) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish()
}



pub fn compute_d_no_min() -> [f64; 128] {
    let mut d = [0.0; 128];

    for i in 0..128 {
        let chr_i = i as u8 as char;
        let ord_i = chr_i as i8;
        let exponent = -(ord_i - 33) as f64 / 10.0;
        d[i] = (10.0_f64).powf(exponent);
    }
    d
}



fn cow_to_string(cow: Cow<'_, [u8]>) -> String {
    String::from_utf8(cow.into_owned()).unwrap_or_else(|e| {
        e.utf8_error().to_string()
        // Handle the error if the conversion fails
    })
}



pub fn get_canonical_kmer_minimizers_hashed(seq: Cow<'_, [u8]>, k_size: usize, w_size: usize, this_minimizers: &mut Vec<Minimizer_hashed>)  {
    //make sure that we have suitable values for k_size and w_size (w_size should be larger)
    let mut w= 0;
    if w_size > k_size{
        w = w_size - k_size;
    }
    //k_size was chosen larger than w_size. To not fail we use every k-mer as minimizer (maybe have an error message?)
    else {
        w = 1;
    }
    //let mut rc_vec=VecDeque::with_capacity(w+1);
    let mut window_kmers: VecDeque<(u64, usize)> = VecDeque::with_capacity(w + 1);
    let mut k_mer_str: &str;
    let mut rc_string;
    let mut forward_hash;
    let mut reverse_hash;
    //let full_seq=cow_to_string(seq.clone());
    //we can only get a minimizer if the sequence is longer than w + k_size - 1 (else we do not even cover one full window)
    if w + k_size < seq.len() + 1{
        for i in 0..w {
            k_mer_str = std::str::from_utf8(&seq[i..i + k_size]).unwrap();
            rc_string = reverse_complement(&k_mer_str);
            //generate the hashes of the kmers
            forward_hash = calculate_hash(&k_mer_str);
            reverse_hash = calculate_hash(&rc_string);
            //we now want to find the canonical minimizer: we only push the smaller k-mer of k_mer_str and rc_String into the window
            if forward_hash <= reverse_hash {
                window_kmers.push_back((forward_hash, i));
            }
            else{
                window_kmers.push_back((reverse_hash, i))
            }
        }
    }
    //println!("kmers in window: {:?}", window_kmers);
    //store the final positional minimizers in a vector
    if !window_kmers.is_empty(){
        // Find the initial minimizer (minimizer of initial window)
        let mut binding= window_kmers.clone();
        let (curr_min, min_pos) = binding.iter().min_by_key(|&(kmer, _)| kmer).unwrap();
        //add the initial minimizer to the vector
        let mut mini = Minimizer_hashed {sequence: *curr_min,position: *min_pos };
        this_minimizers.push(mini.clone());
        //we always store the previous minimizer to compare to the newly found one
        let mut prev_minimizer = mini;
        let mut new_kmer_pos;
        let mut new_kmer_str;
        let mut rc_string;
        let mut forward_hash;
        let mut reverse_hash;
        //iterate further over the sequence and generate the minimizers thereof
        for (i, new_kmer) in seq[w..].windows(k_size).enumerate() {
            new_kmer_pos = i  + w;
            new_kmer_str = std::str::from_utf8(new_kmer).unwrap();
            rc_string = reverse_complement(new_kmer_str);
            // updating  by removing first kmer from window
            window_kmers.pop_front().unwrap();
            forward_hash = calculate_hash(&new_kmer_str);
            reverse_hash = calculate_hash(&rc_string);
            if reverse_hash > forward_hash{
                window_kmers.push_back((forward_hash, new_kmer_pos));
            }
            else {
                window_kmers.push_back((reverse_hash,new_kmer_pos))
            }
            // Find the new minimizer, we need a ds that was cloned from window_kmers to abide ownership rules in rust
            binding = window_kmers.clone();
            let (curr_min, min_pos) = *binding.iter().min_by_key(|&(kmer, _)| kmer).unwrap();
            //make sure that the minimal string is a new minimizer not just the previously found one
            if  min_pos != prev_minimizer.position{ //&& *curr_min != prev_minimizer.1 {
                //add the minimizer into the vector and store the minimizer as previously detected minimizer
                mini = Minimizer_hashed {sequence: curr_min,position: min_pos };
                //println!("minimizer {:?}",mini);
                this_minimizers.push(mini.clone());
                prev_minimizer = mini.clone();
            }
            rc_string.clear();
        }
    }
}




//calculates the average of  a list of f64s and returns it as f64
fn average(numbers: &[f64]) -> f64 {
    numbers.iter().sum::<f64>()/ numbers.len() as f64
}

///Used to detect significant minimizers by checking the read qualities and estimating the overall quality of the area of the read
/// Input: quality_interval: the quality values of the area we want to check
/// Output: significance_indicator: a bool stating whether the minimizer is significant( true: yes, false: no)
///
pub fn is_significant(quality_interval: &[u8], d_no_min:[f64;128],quality_threshold: &f64)->bool{
    let mut significance_indicator= false;
    let mut qualities :Vec<f64> = vec![];
    let mut quality = 1.0;
    let mut index;
    let mut q_value;
    let mut probability_error;
    //for each character in quality string:
    for c in quality_interval {
        index = c.to_ascii_lowercase() as usize;
        //q_value gives the PHRED quality score: i.e. '+' gives us 0.1
        q_value = d_no_min[index];
        //here we get the base call accuracy
        probability_error = 1.0 - q_value;
        qualities.push(probability_error);
        quality *= probability_error
    }

    //TODO: let quality be dependent on length of quality_interval (e.g. 1*E-len)
    if quality > *quality_threshold {
        significance_indicator = true;
    }
    significance_indicator
}


//filter out minimizers for which the quality of the minimizer_impact range is too bad
pub fn filter_minimizers_by_quality(this_minimizers: &Vec<Minimizer_hashed>, fastq_quality:&[u8], k: usize, d_no_min:[f64;128], minimizers_filtered: &mut Vec<Minimizer_hashed>, quality_threshold: &f64) {
    let mut skipped_cter= 0;
    let mut minimizer_range_start;
    let mut significant;
    //println!("Number of minimizers: {}",this_minimizers.len());
    for mini in this_minimizers{
        minimizer_range_start = mini.position;
        let qualitiy_interval = &fastq_quality[minimizer_range_start..minimizer_range_start + k];
        //println!("Quality_interval len {}",qualitiy_interval.len());
        significant = is_significant(qualitiy_interval, d_no_min, quality_threshold);
        //println!("Quality intervallen {}",qualitiy_interval.len());
        if significant{
            minimizers_filtered.push(mini.clone())
        }
        else{
            skipped_cter += 1;
        }
    }
}




///Method used to generate syncmers from reads
/// INPUT:  seq: a string reference to the original read sequence
///         k_size: The size of the k_mer used
///         s_size: The size of s
///         t: The size of parameter t
///OUtput:  syncmers: A vector storing all syncmers (we use the minimizer struct to store them as essentially the same infos)
pub(crate) fn get_kmer_syncmers(seq: &str, k_size: usize, s_size: usize, t: isize) -> Vec<Minimizer> {
    let w = k_size - s_size;

    // get t, the position of s-mer
    // t is chosen to be in the middle of k-mer for the chosen syncmer
    let mut t = t;
    if t < 0 {
        let diff = k_size as isize - s_size as isize;
        t = if diff % 2 == 0 {
            diff / 2
        } else {
            (diff + 1) / 2
        };
        t -= 1;
    }

    let mut syncmers = Vec::new();
    // get list of all s-mers in the first k-mer

    let mut kmer_smers: VecDeque<&str> = (0..=w).map(|i| &seq[i..i + s_size]).collect();
    for i in 0..seq.len() - k_size {
        // add a new syncmer to the list if its smallest s-mer is at place t
        if kmer_smers.iter().position(|&x| x == *kmer_smers.iter().min().unwrap()) == Some(t as usize) {
            syncmers.push(Minimizer {sequence: (seq[ i..i + k_size]).parse().unwrap(), position: i});
        }
        // move the window one step to the right by popping the leftmost
        // s-mer and adding one to the right
        kmer_smers.pop_front();
        kmer_smers.push_back(&seq[i + k_size - s_size + 1..i + k_size + 1]);
    }

    syncmers
}


/*fn syncmers_canonical(seq: &str, k: usize, s: usize, t: usize) -> Vec<Minimizer> {
    let seq_rc = reverse_complement(seq);
    let mut hasher = DefaultHasher::new();
    let mut window_smers_fw: VecDeque<u64> = (0..k - s + 1).map(|i| seq[i..i + s].hash(&mut hasher)).collect();
    let mut window_smers_rc: VecDeque<u64> = (0..k - s + 1).map(|i| seq_rc[i..i + s].hash(&mut hasher)).collect();

    let curr_min_fw = *window_smers_fw.iter().min().unwrap();
    let curr_min_rc = *window_smers_rc.iter().min().unwrap();

    let pos_min_fw = window_smers_fw.iter().position(|&x| x == curr_min_fw).unwrap();
    let pos_min_rc = window_smers_rc.iter().position(|&x| x == curr_min_rc).unwrap();

    let (pos_min, seq_tmp) = if curr_min_fw < curr_min_rc {
        (pos_min_fw, &seq[0..k])
    } else {
        (pos_min_rc, &seq_rc[0..k])
    };

    let mut syncmers:Vec<Minimizer> = vec![];
    if pos_min == t {
        syncmers.push(Minimizer{sequence:seq_tmp.to_string(), position: 0});
    }

    for i in k - s + 1..seq.len() - s {
        seq[i..i + s].hash(&mut hasher);
        let new_smer_fw = hasher.finish();
        seq_rc[i..i + s].hash(&mut hasher);
        let new_smer_rc=hasher.finish();
        // Updating windows
        let _ = window_smers_fw.pop_front();
        window_smers_fw.push_back(new_smer_fw);
        let _ = window_smers_rc.pop_front();
        window_smers_rc.push_back(new_smer_rc);

        let curr_min_fw = *window_smers_fw.iter().min().unwrap();
        let curr_min_rc = *window_smers_rc.iter().min().unwrap();

        let pos_min_fw = window_smers_fw.iter().position(|&x| x == curr_min_fw).unwrap();
        let pos_min_rc = window_smers_rc.iter().position(|&x| x == curr_min_rc).unwrap();

        let (pos_min, seq_tmp) = if curr_min_fw < curr_min_rc {
            (pos_min_fw, &seq[i - (k - s)..i - (k - s) + k])
        } else {
            (pos_min_rc, &seq_rc[i - (k - s)..i - (k - s) + k])
        };

        if pos_min == t {
            let kmer = seq_tmp.to_string();
            syncmers.push(Minimizer {sequence:kmer, position:i - (k - s)});
        }
    }

    syncmers
}*/







#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer_minimizers_0() {
        let input = "ATGCTAGCATGCTAGCATGCTAGC";
        let window_size = 8;
        let k = 3;
        let actual_minimizers = get_kmer_minimizers(input, k, window_size);
        println!("Generated Minimizers: {:?}", actual_minimizers);
        let expected_minimizers = vec![
            Minimizer { sequence: "ATG".to_string(), position: 0 },
            Minimizer { sequence: "AGC".to_string(), position: 5 },
            Minimizer { sequence: "ATG".to_string(), position: 8 },
            Minimizer { sequence: "AGC".to_string(), position: 13 },
            Minimizer { sequence: "ATG".to_string(), position: 16 },
            Minimizer { sequence: "AGC".to_string(), position: 21 },
        ];
        assert_eq!(actual_minimizers, expected_minimizers);
    }
    #[test]
    fn test_kmer_minimizers_1() {
        let input = "CAATTTAAGGCCCGGG";
        let window_size = 10;
        let k = 5;
        let actual_minimizers = get_kmer_minimizers(input, k, window_size);
        println!("Generated Minimizers: {:?}", actual_minimizers);
        let expected_minimizers = vec![
            Minimizer { sequence: "AATTT".to_string(), position: 1 },
            Minimizer { sequence: "AAGGC".to_string(), position: 6 },
            Minimizer { sequence: "AGGCC".to_string(), position: 7 },
        ];
        assert_eq!(actual_minimizers, expected_minimizers);
    }

    #[test]
    fn test_kmer_minimizers_2() {
        let input = "CAAAGTAAGGCCCTCC";
        let window_size = 10;
        let k = 5;
        let actual_minimizers = get_kmer_minimizers(input, k, window_size);
        println!("Generated Minimizers: {:?}", actual_minimizers);
        let expected_minimizers = vec![
            Minimizer { sequence: "AAAGT".to_string(), position: 1 },
            Minimizer { sequence: "AAGGC".to_string(), position: 6 },
            Minimizer { sequence: "AGGCC".to_string(), position: 7 },
        ];
        assert_eq!(actual_minimizers, expected_minimizers);
    }
    #[test]
    fn test_kmer_minimizers_3() {
        let input = "CAATGA";
        let window_size = 10;
        let k = 5;
        let actual_minimizers = get_kmer_minimizers(input, k, window_size);
        println!("Generated Minimizers: {:?}", actual_minimizers);
        let expected_minimizers = vec![
            Minimizer { sequence: "AATGA".to_string(), position: 1 },
        ];
        assert_eq!(actual_minimizers, expected_minimizers);
    }
    #[test]
    fn test_canonical_minis_1(){
        let input ="GGGTAACTTTTCA";
        let window_size=12;
        let k =6;
        let actual_minimizers=get_canonical_kmer_minimizers(input,k,window_size);
        println!("Generated Minimizers: {:?}", actual_minimizers);
        let expected_minimizers = vec![
            Minimizer { sequence: "AAAAGT".to_string(), position: 5 },
        ];
        assert_eq!(actual_minimizers, expected_minimizers);
    }
    #[test]
    fn test_kmer_syncmers(){
        let input ="CATTCAGGAATC";
        let k=5;
        let s=2;
        let t=2;
        let retreived_syncmers=get_kmer_syncmers(input,k,s,t);
        let expected_syncmers=vec![
            Minimizer { sequence: "TCAGG".to_string(), position: 3 },
            Minimizer { sequence: "GGAAT".to_string(),position: 6},

        ];
        assert_eq!(retreived_syncmers, expected_syncmers);
    }
    #[test]
    fn test_average_0(){
        let mut input=vec![];
        input.push(0.5);
        input.push(0.75);
        input.push(0.25);
        let average_res=average(&*input);
        assert_eq!(average_res,0.5);
    }
    #[test]
    fn test_average_1(){
        let mut input=vec![];
        input.push(1.0);
        input.push(2.0);
        input.push(3.0);
        let average_res=average(&*input);
        assert_eq!(average_res,2.0);
    }
}