pub mod ab1;





pub fn add(left: u64, right: u64) -> u64 {
    left + right
}







#[cfg(test)]
mod tests {
    use std::path;

    use crate::ab1::import_ab1;

    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }


    #[test]
    fn test_read_ab1() {
        let fpath = "/data-slow/kangwei-deliver/kangwei-deliver/S22509070002-Epi5A-1.ab1";
        let p = path::Path::new(fpath);
        let records = import_ab1(p).expect("import ab1 error");
        let mut cnt = 0;
        for record in records {
            // println!("{:?}", record);
            cnt += 1;
        }
        println!("{cnt}");



    }
}
