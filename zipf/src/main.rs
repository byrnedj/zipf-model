use std::env;
use rug::{ops::Pow, Float, Assign};


fn compute_normalizer(alpha: Float, m: u32) -> Float {
    let mut total = Float::new(10);
    let mut temp = Float::new(10);
    let mut induction = Float::new(10);

    for i in 1..(m+1) {
        induction.assign(i);
        temp.assign(induction.pow(-1*alpha));
        total.assign(total + temp);
    }
    
    1/total
}


fn populate_subexpressions(m: u32, alpha: Float, norm: Float) -> (Vec<Float>, Vec<Float>) {
    let mut expr: Vec<Float> = Vec::new();
    let mut ln_expr: Vec<Float> = Vec::new();
    let mut induction = Float::new(10);

    for i in 1..(m+1) {
        let mut temp = Float::new(10);
        temp.assign(i);
        expr.push(Float::new(10));
        ln_expr.push(Float::new(10));
        expr[(i-1) as usize].assign(1 - norm*(temp.pow(-1*alpha)));
        ln_expr[(i-1) as usize].assign(-1*expr[(i-1) as usize].ln());
    }

    (expr, ln_expr)
}

fn populate_footprints_derivatives(m: u32, n: u32, expr: Vec<Float>, ln_expr: Vec<Float>) -> (Vec<Float>, Vec<Float>){
    
    let mut fp: Vec<Float> = Vec::new();
    let mut drv: Vec<Float> = Vec::new();

    for i in 0..n {
        fp.push(Float::new(10));
        drv.push(Float::new(10));
        fp[i as usize].clone().assign(0);
        drv[i as usize].clone().assign(0);
        for j in 0..m {
            fp[i as usize].assign(fp[i as usize].clone() + (1 - expr[j as usize].clone().pow(i)));
            drv[i as usize].assign(drv[i as usize].clone() + (ln_expr[j as usize].clone()*(expr[j as usize].clone().pow(i))));
        }
    }

    (fp, drv)
}

fn compute_mrs(m: u32, n: u32, alpha: Float) {
    let norm: Float = compute_normalizer(alpha, m);
    let (expr, ln_expr) : (Vec<Float>, Vec<Float>) = populate_subexpressions(m, alpha, norm);
    let (fp, drv) : (Vec<Float>, Vec<Float>) = populate_footprints_derivatives(m, n, expr, ln_expr);
    println!("parameters: data size = {0}, trace length = {1}, zipf parameter = {2}", m, n, alpha);

    let mut cache_size: u32 = 0;

    for x in 0..(n-1) {
        if fp[x as usize].clone().floor() > cache_size {
            println!("cache size: {0} \n    miss ratio: {1}", fp[x as usize], drv[x as usize]);
            cache_size = fp[x as usize].clone().floor().to_u32_saturating().unwrap();
        }
    }

}


fn main() {
    let args: Vec<String> = env::args().collect();
    let m: u32 = args[1].parse().unwrap();
    let n: u32 = args[2].parse().unwrap();
    let mut alpha = Float::new(10);
    alpha.assign(args[3].parse::<f32>().unwrap());

    compute_mrs(m, n, alpha);

}