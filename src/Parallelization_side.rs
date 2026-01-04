//use bio::io::fastq;
//use bio::io::bed::Record;
//use std::path::Path;
//use rayon::prelude::*;

/*fn split_in_batches(n_threads: usize, score_vec: Vec<(i32,usize)>, reads_split: &mut Vec<Vec<i32>>) {
    let each_len = score_vec.len() / n_threads + if score_vec.len() % n_threads == 0 {0} else {1};
    //let mut out = vec![Vec::with_capacity(each_len); n_threads];
    //let mut reader = fastq::Reader::from_file(Path::new(filename)).expect("We expect the file to exist");
    for read_info in score_vec.into_iter() {
        //for(i,d) in data.iter().copied().enumerate(){
       reads_split[i % n_threads].push(d);
   }
    //}

}*/
