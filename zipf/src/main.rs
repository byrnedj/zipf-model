use std::env;
use math::round;
use rug::{ops::Pow, Float, Assign};


fn compute_normalizer(alpha: &Float, m: &u32) -> Float {
    let mut total = Float::new(100);
    let mut temp = Float::new(100);
    let mut induction = Float::new(100);

    for i in 1..(m+1) {
        induction.assign(i);
        temp.assign((&induction).pow(-1*(alpha.clone())));
        total.assign(total.clone() + temp.clone());
    }
    let norm = 1.0/total;
    
    norm
}


fn populate_subexpressions(m: &u32, alpha: &Float, norm: &Float) -> (Vec<Float>, Vec<Float>, Vec<Float>) {
    let mut expr: Vec<Float> = Vec::new();
    let mut ln_expr: Vec<Float> = Vec::new();
    let mut wfp_expr: Vec<Float> = Vec::new();


    //------------------change rw ratio here---------------------
    let rw = 0.9;
    //-----------------------------------------------------------

    for i in 1..(m+1) {
        let mut temp = Float::new(100);
        temp.assign(i);
        expr.push(Float::new(100));
        ln_expr.push(Float::new(100));
        wfp_expr.push(Float::new(100));
        expr[(i-1) as usize].assign(1 - norm*(temp.clone().pow(-1*(alpha.clone()))));
        //println!("{}", &expr[(i-1) as usize]);
        ln_expr[(i-1) as usize].assign(-1*expr[(i-1) as usize].clone().ln());
        wfp_expr[(i-1) as usize].assign((1.0-rw).pow(norm*(temp.clone().pow(-1*(alpha.clone())))));
    }

    (expr, ln_expr, wfp_expr)
}

fn populate_footprints_derivatives(m: &u32, n: &u32, expr: Vec<Float>, ln_expr: Vec<Float>, wfp_expr: Vec<Float>) -> (Vec<Float>, Vec<Float>, Vec<Float>){
    
    let mut fp: Vec<Float> = Vec::new();
    let mut drv: Vec<Float> = Vec::new();
    let mut wfp: Vec<Float> = Vec::new();
    let mut rdfp: Vec<Float> = Vec::new();

    for i in 0..*n {
        fp.push(Float::new(100));
        drv.push(Float::new(100));
        wfp.push(Float::new(100));
        rdfp.push(Float::new(100));

        fp[i as usize].assign(0);
        drv[i as usize].assign(0);

        for j in 0..*m {
            let t1 = fp[i as usize].clone();
            let t2 = drv[i as usize].clone();
            fp[i as usize].assign(t1 + (1 - expr[j as usize].clone().pow(i)));
            drv[i as usize].assign(t2 + (ln_expr[j as usize].clone()*(expr[j as usize].clone().pow(i))));
            
        }
        
        wfp[i as usize].assign(fp[i as usize].clone());

        for k in 0..*m{
            let t1 = wfp_expr[k as usize].clone();
            let t2 = 1 - expr[k as usize].clone();
            let t3 = wfp[i as usize].clone();

            wfp[i as usize].assign(t3 - (t1.pow(i)*t2));
        }

        let t1 = fp[i as usize].clone();
        let t2 = wfp[i as usize].clone();

        rdfp[i as usize].assign(t2/t1);
    }

    (fp, drv, rdfp)
}

fn compute_mrs(m: u32, n: u32, alpha: Float) {
    let norm: Float = compute_normalizer(&alpha, &m);
    let (expr, ln_expr, wfp_expr) : (Vec<Float>, Vec<Float>, Vec<Float>) = populate_subexpressions(&m, &alpha, &norm);
    let (fp, drv, rdfp) : (Vec<Float>, Vec<Float>, Vec<Float>) = populate_footprints_derivatives(&m, &n, expr, ln_expr, wfp_expr);
    //println!("parameters: data size = {0}, trace length = {1}, zipf parameter = {2}", &m, &n, &alpha);

    let mut cache_size: u32 = 0;
    println!("size,mr,wbr");
    for x in 0..(n-1) {
        if fp[x as usize].clone().floor() > cache_size {
            let size = round::floor(fp[x as usize].to_f64(),1);
            if size.fract() == 0.0 {
                let mr = round::floor(drv[x as usize].to_f64(),3);
                let wbr = round::floor(rdfp[x as usize].to_f64(),3) * mr;
                println!("{0},{1:.3},{2:.3}",size,mr,wbr);
            }
            cache_size = fp[x as usize].clone().floor().to_u32_saturating().unwrap();
        }
    }

}


fn main() {
    let args: Vec<String> = env::args().collect();
    let m: u32 = args[1].parse().unwrap();
    let n: u32 = args[2].parse().unwrap();
    let mut alpha = Float::new(100);
    alpha.assign(args[3].parse::<f32>().unwrap());

    compute_mrs(m, n, alpha);

}
