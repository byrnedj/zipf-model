use std::env;
use rug::{Integer, Float}


fn compute_normalizer(alpha: Float, m: u32) -> Float {
    let mut total = Float::new();
    let mut temp = Float::new();
    let mut induction = Float::new();

    for i in 1..(m+1) {
        induction.assign(i);
        temp.assign(induction.pow(-1*alpha))
        total.assign(total + temp)
    }
    
    1/total
}


fn populate_subexpressions(m: u32, alpha: Float, norm: Float) -> ([Float; m], [Float; m]) {
    let mut expr: [Float, m];
    let mut ln_expr: [Float, m];
    let mut induction = Float::new();

    for i in 1..(m+1) {
        expr[i-1].assign(1 - norm*(i.pow(-1*alpha)));
        ln_expr[i-1].assign(-1*expr[i-1].ln())
    }

    (expr, ln_expr)
}

fn populate_footprints_derivatives(m: u32, n: u32, expr: [Float; m], ln_expr: [Float; m]) -> ([Float; n], [Float; n]){
    
    let mut fp: [Float, n];
    let mut drv: [Float, n];

    for i in 0..n {
        fp[i].assign(0);
        drv[i].assign(0);
        for j in 0..m {
            fp[i].assign(fp[i] + (1 - expr[j].pow(i)))
            drv[i].assign(drv[i] + (ln_expr[j]*(expr[j].pow(i))))
        }
    }

    (fp, drv)
}

fn compute_mrs(m: u32, n: u32, alpha: f32) {
    let norm: f32 = compute_normalizer(alpha, m);
    let (expr, ln_expr) : ([Float; m], [Float; m]) = populate_subexpressions(m, alpha, norm);
    let (fp, drv) : ([Float; n], [Float; n]) = populate_footprints_derivatives(m, n, expr, ln_expr);
    println!("parameters: data size = {0}, trace length = {1}, zipf parameter = {2}", m, n, alpha);

    let mut cache_size: u32 = 0;

    for x in 0..(n-1):
        if fp[x].floor() > cache_size {
            println!("cache size: {0} \n    miss ratio: {1}", fp[x], drv[x]);
            cache_size = math.floor(fp[x]);
        }

}


fn main() {
    let args: Vec<String> = env::args().collect();
    let m: u32 = &args[1].parse().unwrap();
    let n: u32 = &args[2].parse().unwrap();
    let alpha = Float::new();
    alpha.assign(&args[3].parse().unwrap());

    compute_mrs(m, n, alpha);

}