use core::f64;
use std::fs::File;
use std::io::Write;
use num::pow::Pow;
use num::Complex;
use std::f64::consts::PI;
use rayon::prelude::*;

fn main() {
    let n_points: usize = 16000;
    let dt: f64 = 0.002;
    let n_steps: usize = 2000;
    let epsilon: f64 = 0.05;  

    let mut z = init_z(n_points,1.0);
    let gamma = init_gamma(n_points,1.0);

    let mut file = std::fs::File::create("data.csv").unwrap();

    let mut nudge: Vec<Complex<f64>> = Vec::new();
    nudge.resize(n_points, Complex::new(0.0,0.0));

    for i in 0..n_steps
    {
        println!("{}",i);
        //if i % 10 == 0 {
            print(&z,&nudge, n_points,i, &mut file);
        //}

        z = rk4(&z,&gamma,n_points,dt,epsilon);

        for j in 0..n_points {
            if z[j].re > 1.0 {
                z[j] -= 1.0;
                nudge[j] += 1.0;
            } else if z[j].re < 0.0 {
                z[j] += 1.0;
                nudge[j] -= 1.0;
            }
        }
        
    }

}

fn print(z: &[Complex<f64>], nudge: &[Complex<f64>], n: usize, t: usize, file: &mut File) 
{
    for i in 0..n {
        let _ = writeln!(file,"{},{},{},{}", t, i, z[i].re + nudge[i].re, z[i].im);
    }
}

fn init_z(n: usize, k: f64) -> Vec<Complex<f64>>
{
    let d_xi = (n as f64).recip();
    let mut xi = 0.0;

    let mut answer: Vec<Complex<f64>> = Vec::new();

    for _i in 1..=n/2 
    {
        xi = xi + 2.0 * d_xi;
         answer.push(Complex::new(xi,0.2 * (2.0 * k * PI * xi).sin()));
        //answer.push(Complex::new(xi + 0.1 * (2.0 * k * PI * xi).sin(),0.2 * (2.0 * k * PI * xi).sin()));
        /*answer.push(
            Complex::new(
                xi + 0.05 * (2.0 * 2.0 * k * PI * xi).sin(),
                0.05 * (2.0 * 2.0 * k * PI * xi).sin() * if xi > 0.5 { -1.0 } else { 1.0 }
            )
        );*/

        //answer.push(Complex::new(xi, -0.2));
    }

    let i1 : Complex<f64> = Complex::new(0.0, -0.5); // original -0.2; thin -0.1; thick -0.5

    let temp : Vec<Complex<f64>> = answer.clone().iter().map(|z| z + i1).collect();
    temp.iter().for_each(|z| answer.push(*z));
    answer
}

fn init_gamma(n: usize, k: f64) -> Vec<f64>
{
    let d_xi = (n as f64).recip();
    let mut xi = 0.0;

    let mut answer: Vec<f64> = Vec::new();

    for _i in 1..=n 
    {
        xi = xi + d_xi;
        //answer.push((2.0 * k * PI * xi).sin());
        //answer.push((2.0 * 2.0 * k * PI * xi).sin() * if xi > 0.5 { -1.0 } else { 1.0 });
        answer.push(1.0);
        answer.push(1.0);
    }
    answer
}

fn rk4(z: &[Complex<f64>], gamma: &[f64], n: usize, h: f64, epsilon: f64) -> Vec<Complex<f64>>
{
    let k1 : Vec<Complex<f64>> = f_z(z,gamma,n,epsilon).iter().map(|u| u * h).collect();
    let k0_k12: Vec<Complex<f64>> = z.iter().zip(k1.iter()).map(|(a,b)| a + b / 2.0).collect();
    let k2 : Vec<Complex<f64>> = f_z(&k0_k12, gamma,n,epsilon).iter().map(|u| u * h).collect();
    let k0_k22 : Vec<Complex<f64>> = z.iter().zip(k2.iter()).map(|(a,b)| a + b / 2.0).collect();
    let k3 : Vec<Complex<f64>> = f_z(&k0_k22,gamma,n,epsilon).iter().map(|u| u * h).collect();
    let k0_k3 : Vec<Complex<f64>> = z.iter().zip(k3.iter()).map(|(a,b)| a + b).collect();
    let k4 : Vec<Complex<f64>> = f_z(&k0_k3,gamma,n,epsilon).iter().map(|u| u * h).collect();

    k1.iter().zip(k2.iter()).map(|(a,b)| a + 2.0 * b)
        .zip(k3.iter()).map(|(a,b)| a + 2.0 * b)
        .zip(k4.iter()).map(|(a,b)| (a + b) / 6.0) 
        .zip(z.iter()).map(|(a,b)| a + b).collect()
}

fn f_z(z: &[Complex<f64>], gamma: &[f64], n: usize, epsilon: f64) -> Vec<Complex<f64>>
{
    let mut answer: Vec<Complex<f64>> = Vec::new();

    (0..n)
        .into_par_iter()
        .map(|i|
        {
            let mut accum: Complex<f64> = Complex::new(0.0,0.0);
            for j in 0..n
            {
                if i == j { 
                   continue;
                }
                
                let dz = z[i] - z[j];

                let temp_real = (-2.0 * PI * dz.im).sinh();
                let temp_im = (2.0 * PI * dz.re).sin();
                let temp_denom = (dz.im * 2.0 * PI).cosh() - (dz.re * 2.0 * PI).cos() + epsilon.pow(2);

                accum += Complex::new(temp_real / temp_denom, temp_im / temp_denom) * gamma[j];
            }
            accum /= 2.0 * (n as f64);
            accum
        }
    )
    .collect_into_vec(&mut answer);
    answer
}