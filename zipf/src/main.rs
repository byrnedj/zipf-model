use math::round;
use rug::{ops::Pow, Float, Assign};
extern crate rgsl;
extern crate clap;
use clap::{Arg, App};

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
    let mut binom_expr: Vec<Float> = Vec::new();

    for i in 1..(m+1) {
        let mut temp = Float::new(100);
        temp.assign(i);
        expr.push(Float::new(100));
        ln_expr.push(Float::new(100));
        binom_expr.push(Float::new(100));
        expr[(i-1) as usize].assign(1 - norm*(temp.clone().pow(-1*(alpha.clone()))));
        ln_expr[(i-1) as usize].assign(-1*expr[(i-1) as usize].clone().ln());
        binom_expr[(i-1) as usize].assign(norm*(temp.clone().pow(-1*(alpha.clone()))));
    }

    (expr, ln_expr, binom_expr)
}

fn populate_footprints_derivatives(m: &u32, n: &u32, alpha: &Float, norm:&Float, r: &f32, expr: Vec<Float>, ln_expr: Vec<Float>, binom_expr: Vec<Float>) -> (Vec<Float>, Vec<Float>, Vec<Float>){
    
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

        wfp[i as usize].assign(0);
        for k in 0..*m{
            let t1 = binom_expr[k as usize].clone()*i;
            let t3 = wfp[i as usize].clone();
 
            wfp[i as usize].assign(t3 + (1 - ((1.0-r).pow(t1))));
        }
        

        let t1 = fp[i as usize].clone();
        let t2 = wfp[i as usize].clone();

        rdfp[i as usize].assign(t2/t1);
    }

    (fp, drv, rdfp)
}

fn populate_approx(m: &u32, n: &u32, alpha: &Float, norm: &Float, r: &f32) -> (Vec<Float>, Vec<Float>, Vec<Float>){
    let mut fp_approx: Vec<Float> = Vec::new();
    let mut drv_approx: Vec<Float> = Vec::new();
    let mut wfp_approx: Vec<Float> = Vec::new();
    let mut rdfp_approx: Vec<Float> = Vec::new();

    let inva: f64 = 1.0/(alpha.to_f64());

    for i in 0..*n{
        fp_approx.push(Float::new(100));
        drv_approx.push(Float::new(100));
        wfp_approx.push(Float::new(100));
        rdfp_approx.push(Float::new(100));

        //TODO: fp approx, drv approx
        let t: f64 = norm.to_f64() * (i as f64)/(*m as f64).pow(alpha.to_f64());
        let fp_temp: f64 = (*m as f64)*(1.0 - (inva*t.pow(inva)*rgsl::gamma_beta::incomplete_gamma::gamma_inc(-1.0*inva, t)));
        let drv_temp: f64 = ((i as f64*norm.to_f64()).pow(inva)*rgsl::gamma_beta::incomplete_gamma::gamma_inc(1.0 - inva, t))/(alpha.to_f64()*i as f64);
        fp_approx[i as usize].assign(fp_temp);
        drv_approx[i as usize].assign(drv_temp);
        

        //wfp approximation using definite integral
        wfp_approx[i as usize].assign(m);

        let mut norm_temp = Float::new(100);
        let mut temp1 = Float::new(100);
        let mut tempM = Float::new(100);
        let mut exp1 = Float::new(100);
        let mut expM = Float::new(100);
        let mReciprocal: f64 = 1.0/(*m as f64);
        let mut expMtemp = Float::new(100);
        let n_func: f64 = 1.0 + (1.0/alpha.to_f64());

        expMtemp.assign(mReciprocal.pow(alpha));
        norm_temp.assign(norm*(-1.0*(i as f64)));
        temp1.assign(norm_temp.clone()*((1.0-r).ln()));
        tempM.assign(norm_temp.clone()*expMtemp*((1.0-r).ln()));

        exp1.assign((temp1.clone().to_f64().pow(n_func - 1.0))*rgsl::gamma_beta::incomplete_gamma::gamma_inc(1.0 - n_func, temp1.to_f64()));
        //exp1.assign(rgsl::gamma_beta::incomplete_gamma::gamma_inc(1.0 - n_func, temp1.to_f64()));
        expM.assign(((*m as f64)*tempM.clone().pow(n_func - 1.0)) * rgsl::gamma_beta::incomplete_gamma::gamma_inc(1.0 - n_func, tempM.to_f64()));

        wfp_approx[i as usize].assign(m - 1 - (expM.clone()- exp1.clone())/alpha);

        let t1 = fp_approx[i as usize].clone();
        let t2 = wfp_approx[i as usize].clone();

        rdfp_approx[i as usize].assign(t2/t1);
    }

    (fp_approx, drv_approx, rdfp_approx)
}

fn compute_approx(m: u32, n: u32, alpha: Float, r: f32, output_int: u32) {
    let norm: Float = compute_normalizer(&alpha, &m);
    //let (expr, ln_expr, binom_expr) : (Vec<Float>, Vec<Float>, Vec<Float>) = populate_subexpressions(&m, &alpha, &norm);
    let (fp_approx, drv_approx, rdfp_approx) : (Vec<Float>, Vec<Float>, Vec<Float>) = populate_approx(&m, &n, &alpha, &norm, &r);

    let mut cache_size: u32 = 0;
    println!("size,mr,wbr");
    for x in 1..(n-1) {
        let fp = fp_approx[x as usize].clone().floor();
        if  fp > cache_size && fp > output_int {
            let fsize = output_int as f64;
            let size = round::floor(fp_approx[x as usize].to_f64(), fsize.log10() as u8) as u32;
            let isize = (size/output_int)*output_int;
            let mr = round::floor(drv_approx[x as usize].to_f64(), 3);
            let wbr = round::floor(rdfp_approx[x as usize].to_f64(), 3) * mr;
            println!("{0},{1:.3},{2:.3}",isize,mr,wbr);
            cache_size = fp_approx[x as usize].clone().floor().to_u32_saturating().unwrap() + output_int;
        }
    }
}

fn compute(m: u32, n: u32, alpha: Float, r: f32) {
    let norm: Float = compute_normalizer(&alpha, &m);
    let (expr, ln_expr, binom_expr) : (Vec<Float>, Vec<Float>, Vec<Float>) = populate_subexpressions(&m, &alpha, &norm);
    let (fp, drv, rdfp) : (Vec<Float>, Vec<Float>, Vec<Float>) = populate_footprints_derivatives(&m, &n, &alpha, &norm, &r, expr, ln_expr, binom_expr);

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
    //clap definition
    let matches = App::new("zipfmodel")
    .author("Wesley Smith <wsmith6@ur.rochester.edu>")
    .arg(Arg::with_name("type")
             .short("t")
             .long("type")
             .takes_value(true)
             .required(true)
             .help("exact or approximate math. vals = exact, approx"))
    .arg(Arg::with_name("datasize")
             .short("m")
             .long("datasize")
             .takes_value(true)
             .required(true)
             .help("Number of distinct elements"))
    .arg(Arg::with_name("tracelength")
             .short("n")
             .long("tracelength")
             .takes_value(true)
             .required(true)
             .help("Maximum window size"))
    .arg(Arg::with_name("alpha")
             .short("a")
             .long("alpha")
             .takes_value(true)
             .required(true)
             .help("Zipf parameter"))
    .arg(Arg::with_name("write_ratio")
             .short("r")
             .long("write_ratio")
             .takes_value(true)
             .required(true)
             .help("percentage of accesses that are writes"))
    .arg(Arg::with_name("output_int")
             .short("i")
             .long("output_int")
             .takes_value(true)
             .required(true)
             .help("MRC/WBR interval"))
    .get_matches();
    
    //parse clap clargs
    let ex_type = matches.value_of("type").unwrap();
    let mut alpha = Float::new(100);
    let a_str: &str = matches.value_of("alpha").unwrap();
    let a = a_str.parse::<f32>().unwrap();
    alpha.assign(a);
    let m_str: &str = matches.value_of("datasize").unwrap();
    let m: u32 = m_str.parse::<u32>().unwrap();
    let n_str: &str = matches.value_of("tracelength").unwrap();
    let n = n_str.parse::<u32>().unwrap();
    let r_str: &str = matches.value_of("write_ratio").unwrap();
    let r = r_str.parse::<f32>().unwrap();
    let output_int_str: &str = matches.value_of("output_int").unwrap();
    let output_int = output_int_str.parse::<u32>().unwrap();

    //error handler off
    rgsl::error::set_error_handler_off();

    //execute based on type
    if ex_type.eq("exact"){
        compute(m, n, alpha, r);
    }
    else if ex_type.eq("approx"){
        compute_approx(m, n, alpha, r, output_int);
    }
    else{
        println!("Invalid execution type. cargo run -h for help");
    }

}
