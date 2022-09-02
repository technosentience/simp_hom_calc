#![feature(test)]
#![feature(is_some_with)]

mod homology;
mod parser;
mod simplex;
mod smith;
mod test;

use homology::homology_groups;
use parser::parse_complex;
use simplex::SimplicialComplex;
use std::io;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let input = io::stdin();
    let mut buf = String::new();

    while parse_complex(&buf).is_err_and(nom::Err::is_incomplete) {
        let mut line = String::new();
        input.read_line(&mut line)?;
        buf += &line;
    }

    let (_, simplices) =
        nom::Finish::finish(parse_complex(buf.as_str())).map_err(|s| s.to_string())?;
    let simplices = SimplicialComplex::new(simplices);
    let homs = homology_groups(simplices);
    for (dim, hom) in homs.iter().enumerate() {
        println!("H_{}: {}", dim, hom)
    }

    return Ok(());
}
