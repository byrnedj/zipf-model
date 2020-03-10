use std::env;
//TODO: GMP crate

fn compute_normalizer(alpha: f32, m: u32) -> f32 {
    //TODO
}

fn populate_derivatives(m: u32, n: u32) {
    //TODO
}

fn populate_subexpressions(m: u32, alpha: f32, norm: f32) {
    //TODO
}

fn populate_footprint(m: u32, n: u32) {
    //TODO
}

fn compute_mrs(m: u32, n: u32, alpha: f32) {
    let norm: f32 = compute_normalizer(alpha, m);
    populate_subexpressions(m, alpha, norm);
    populate_footprint(m, n);
    populate_derivatives(m, n);
    println!("parameters: data size = {0}, trace length = {1}, zipf parameter = {2}", m, n, alpha);

    let cache_size: u32 = 0;

    //TODO
}


fn main() {
    let args: Vec<String> = env::args().collect();
    let m: u32 = &args[1].parse().unwrap();
    let n: u32 = &args[2].parse().unwrap();
    let alpha: f32 = &args[3].parse().unwrap();
    compute_mrs(m, n, alpha);


}