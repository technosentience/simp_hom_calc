use nom::{
    character::streaming::char, character::streaming::multispace0, character::streaming::u64,
    combinator::map, error::ParseError, multi::separated_list1, sequence::delimited, IResult,
};

// Skip whitespaces
fn ws<'a, F: 'a, O, E: ParseError<&'a str>>(
    inner: F,
) -> impl FnMut(&'a str) -> IResult<&'a str, O, E>
where
    F: FnMut(&'a str) -> IResult<&'a str, O, E>,
{
    delimited(multispace0, inner, multispace0)
}

fn parse_list<'a, F: 'a, O, E: ParseError<&'a str> + 'a>(
    inner: F,
) -> impl FnMut(&'a str) -> IResult<&'a str, Vec<O>, E>
where
    F: FnMut(&'a str) -> IResult<&'a str, O, E>,
{
    delimited(
        ws(char('[')),
        separated_list1(ws(char(',')), inner),
        char(']'),
    )
}

fn parse_simplex(input: &str) -> IResult<&str, Vec<usize>> {
    parse_list(map(ws(u64), |x| x as usize))(input)
}

pub fn parse_complex(input: &str) -> IResult<&str, Vec<Vec<usize>>> {
    parse_list(ws(parse_simplex))(input)
}
