//! Writes the model `reduce_with_presolve` hands the search as LP, to check it elsewhere.
use ripsolve::{Problem, reduce};

fn main() {
    let mut args = std::env::args().skip(1);
    let (input, output) = (args.next().unwrap(), args.next().unwrap());
    let original = Problem::from_file(std::path::Path::new(&input)).unwrap();
    match reduce::reduce_with_presolve(&original, 1e-9) {
        Ok(Some((smaller, map, stats))) => {
            println!("{} cols kept, offset {}, {stats:?}", map.columns().len(), smaller.obj_offset);
            std::fs::write(output, smaller.to_lp()).unwrap();
        }
        other => println!("{:?}", other.map(|o| o.is_some())),
    }
}
