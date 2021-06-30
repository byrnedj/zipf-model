use math::round;
use rug::{ops::Pow, Float, Assign};
extern crate rgsl;
extern crate clap;
use clap::{Arg, App};
use std::fs::File;
use std::io::Write;

fn compute_normalizer(alpha: &f64, m: &u32) -> f64 {
    let mut total: f64 = 0.0;
    let mut temp: f64 = 0.0;
    let mut induction: f64 = 0.0;

    for i in 1..(m+1) {
        induction = i as f64;
        temp = (&induction).pow(-1.0*(alpha.clone()));
        total = total + temp;
    }
    let norm = 1.0/total;
    
    norm
}


fn populate_subexpressions(m: &u32, alpha: &f64, norm: &f64) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut expr: Vec<f64> = Vec::new();
    let mut ln_expr: Vec<f64> = Vec::new();
    let mut binom_expr: Vec<f64> = Vec::new();

    for i in 1..(m+1) {
        let mut temp: f64 = 0.0;
        temp = i as f64;
        expr.push(0.0);
        ln_expr.push(0.0);
        binom_expr.push(0.0);
        expr[(i-1) as usize] = 1.0 - norm*(temp.clone().pow(-1.0*(alpha.clone())));
        ln_expr[(i-1) as usize] = -1.0*expr[(i-1) as usize].clone().ln();
        binom_expr[(i-1) as usize] = norm*(temp.clone().pow(-1.0*(alpha.clone())));
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

fn populate_approx(m: &u32, n: &u32, alpha: &f64, norm: &f64, r: &f32, c: u32) -> (f64, f64, f64){
    let mut fp_approx: f64 = 0.0;
    let mut drv_approx: f64 = 0.0;
    let mut wfp_approx: f64 = 0.0;
    let mut rdfp_approx: f64 = 0.0;

    let inva: f64 = 1.0/alpha;
    let mut searchup = false;
    let mut searchdown = false;
    let mut win_length: f64 = 1.0;
    let mut lower_bound: f64 = 1.0;
    let mut t_final: f64 = 0.0;

    while !searchup{
        let t: f64 = norm * (win_length as f64)/(*m as f64).pow(alpha);
        let fp_temp: f64 = (*m as f64)*(1.0 - (inva*t.pow(inva)*rgsl::gamma_beta::incomplete_gamma::gamma_inc(-1.0*inva, t)));
        if fp_temp > (c as f64){
            lower_bound = win_length/2.0;
            searchup = true;
        }
        else{
            win_length = win_length * 2.0;
        }
    }

    while !searchdown{
        let t: f64 = norm * ((win_length + lower_bound)/2.0)/(*m as f64).pow(alpha);
        let t2: f64 = norm * (((win_length + lower_bound)/2.0) - 1.0)/(*m as f64).pow(alpha);
        let fp_temp: f64 = (*m as f64)*(1.0 - (inva*t.pow(inva)*rgsl::gamma_beta::incomplete_gamma::gamma_inc(-1.0*inva, t)));
        let fp_temp2: f64 = (*m as f64)*(1.0 - (inva*t.pow(inva)*rgsl::gamma_beta::incomplete_gamma::gamma_inc(-1.0*inva, t2)));

        if fp_temp > (c as f64) && fp_temp2 < (c as f64){
            searchdown = true;
            t_final = t;
            win_length = (win_length + lower_bound)/2.0;
        }
        else if fp_temp > (c as f64){
            win_length = (win_length + lower_bound)/2.0;
        }
        else{
            lower_bound = (win_length + lower_bound)/2.0;
        }
    }

    //TODO: fp approx, drv approx
    fp_approx = (*m as f64)*(1.0 - (inva*t_final.pow(inva)*rgsl::gamma_beta::incomplete_gamma::gamma_inc(-1.0*inva, t_final)));
    drv_approx = ((win_length as f64*norm).pow(inva)*rgsl::gamma_beta::incomplete_gamma::gamma_inc(1.0 - inva, t_final))/(alpha*win_length as f64);


    //wfp approximation using definite integral
    let mut norm_temp: f64 = 0.0;
    let mut temp1: f64 = 0.0;
    let mut tempM: f64 = 0.0;
    let mut exp1: f64 = 0.0;
    let mut expM: f64 = 0.0;
    let mReciprocal: f64 = 1.0/(*m as f64);
    let mut expMtemp: f64 = 0.0;
    let n_func: f64 = 1.0 + (1.0/alpha);

    expMtemp = mReciprocal.pow(alpha);
    norm_temp = norm*(-1.0*(win_length as f64));
    temp1 = norm_temp.clone()*((1.0-(*r as f64)).ln());
    tempM = norm_temp.clone()*expMtemp*((1.0-(*r as f64)).ln());

    exp1 = (temp1.clone().pow(n_func - 1.0))*rgsl::gamma_beta::incomplete_gamma::gamma_inc(1.0 - n_func, temp1);
    expM =((*m as f64)*tempM.clone().pow(n_func - 1.0)) * rgsl::gamma_beta::incomplete_gamma::gamma_inc(1.0 - n_func, tempM);

    wfp_approx = (*m as f64) - 1.0 - (expM.clone()- exp1.clone())/alpha;

    rdfp_approx = wfp_approx/fp_approx.clone();
    
    (fp_approx, drv_approx, rdfp_approx)
}

fn compute_approx(m: u32, n: u32, alpha: f64, r: f32, output_int: u32, use_stdout: u32, ofilename: &str, c: u32, norm: f64) {
    let (fp_approx, drv_approx, rdfp_approx) : (f64, f64, f64) = populate_approx(&m, &n, &alpha, &norm, &r, c);
    let mut of;
    if use_stdout == 0 {
        of = File::create(ofilename).expect("ERROR!");
    } else {
        of = File::create("/dev/null").unwrap();
    }

    //let mut cache_size: u32 = 0;
    //if use_stdout == 0 {
    //    let line = format!("size,mr,wbr\n");
    //    of.write(line.as_bytes()).expect("write ERROR");
    //} else {
    //    println!("size,mr,wbr");
    //}
    println!("size,mr,wbr");
    
    let fsize = output_int as f64;
    let size = round::floor(fp_approx, fsize.log10() as u8) as u32;
    let isize = (size/output_int)*output_int;
    let mr = round::floor(drv_approx, 3);
    let wbr = round::floor(rdfp_approx, 3) * mr;
    //if use_stdout == 0 {
    //    let line = format!("{0},{1:.3},{2:.3}\n",isize,mr,wbr);
    //    of.write(line.as_bytes()).expect("WRITE_ERROR");
    //} else {
    //    println!("{0},{1:.3},{2:.3}",isize,mr,wbr);
    //}
    //cache_size = fp_approx.clone().floor().to_u32_saturating().unwrap() + output_int;
    println!("{0},{1:.3},{2:.3}",c,mr,wbr);
        
    
}

fn compute(m: u32, n: u32, alpha: f64, r: f32) {
//    let norm: Float = compute_normalizer(&alpha, &m);
//    let (expr, ln_expr, binom_expr) : (Vec<Float>, Vec<Float>, Vec<Float>) = populate_subexpressions(&m, &alpha, &norm);
//    let (fp, drv, rdfp) : (Vec<Float>, Vec<Float>, Vec<Float>) = populate_footprints_derivatives(&m, &n, &alpha, &norm, &r, expr, ln_expr, binom_expr);

//    let mut cache_size: u32 = 0;
//    for x in 0..(n-1) {
//        if fp[x as usize].clone().floor() > cache_size {
//            let size = round::floor(fp[x as usize].to_f64(),1);
//            if size.fract() == 0.0 {
//                let mr = round::floor(drv[x as usize].to_f64(),3);
//                let wbr = round::floor(rdfp[x as usize].to_f64(),3) * mr;
//                println!("{0},{1:.3},{2:.3}",size,mr,wbr);
//            }
//            cache_size = fp[x as usize].clone().floor().to_u32_saturating().unwrap();
//        }
//    }
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
    .arg(Arg::with_name("output_file")
             .short("o")
             .long("output_file")
             .takes_value(true)
             .required(false)
             .help("outputs MRC/WBR to file"))
    .get_matches();

    let mut use_stdout = 0;
    let mut ofilename = "";

    //parse clap clargs
    let ex_type = matches.value_of("type").unwrap();
    let mut alpha: f64 = 0.0;
    let a_str: &str = matches.value_of("alpha").unwrap();
    let a = a_str.parse::<f64>().unwrap();
    alpha = a;
    let m_str: &str = matches.value_of("datasize").unwrap();
    let m: u32 = m_str.parse::<u32>().unwrap();
    let n_str: &str = matches.value_of("tracelength").unwrap();
    let n = n_str.parse::<u32>().unwrap();
    let r_str: &str = matches.value_of("write_ratio").unwrap();
    let r = r_str.parse::<f32>().unwrap();
    let output_int_str: &str = matches.value_of("output_int").unwrap();
    let output_int = output_int_str.parse::<u32>().unwrap();
    
    let output_file = matches.value_of("output_file");
    match output_file {
        None => use_stdout = 1,
        Some(s) => ofilename = s
    }

    let mut cache_sizes: Vec<u32> = Vec::new();
    for x in 1..2u32 {
        println!("{0}",x*10000);
        cache_sizes.push(x*10000);
    }
 
    //let cache_sizes: [u32; 4] = [1000000,1500000, 2000000];

    //error handler off
    rgsl::error::set_error_handler_off();

    
    //execute based on type
    if ex_type.eq("exact"){
        compute(m, n, alpha, r);
    }

    else if ex_type.eq("approx"){
        let norm: f64 = compute_normalizer(&alpha, &m);
        println!("size,mr,wbr");
        for c in 0..cache_sizes.len(){
            compute_approx(m, n, alpha.clone(), r, output_int, use_stdout, ofilename, cache_sizes[c], norm.clone());
        }
        
    }
    else{
        println!("Invalid execution type. cargo run -h for help");
    }

}
