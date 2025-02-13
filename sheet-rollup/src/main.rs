use std::{f64, fs::File};
use std::io::Write;
use std::thread;

fn main() {
    let n_points: usize = 800;
    let dt: f64 = 0.0005;
    let n_steps: usize = 1000;
    let epsilon: f64 = 0.05;  

    let mut x = init_x(n_points);
    let mut y = init_y(n_points);

    let mut file = std::fs::File::create("data.csv").unwrap();

    for i in 0..n_steps
    {
        println!("{}",i);
        //if i % 10 == 0 {
            print(&x,&y,n_points,i, &mut file);
        //}

        let (new_x, new_y) = thread::scope(|s| {

            let join1 = s.spawn(|| rk4_x(&x,&y,n_points,dt,epsilon));
            let join2 = s.spawn(|| rk4_y(&x,&y,n_points,dt,epsilon));

            (join1.join(), join2.join())

        });

        x = new_x.unwrap();
        y = new_y.unwrap();
        
    }

}

fn print(x: &[f64], y: &[f64], n: usize, t: usize, file: &mut File) 
{
    for i in 0..n {
            let _ = writeln!(file,"{},{},{},{}", t, i, x[i], y[i]);
    }
}

fn step_x(x: &Vec<f64>, y: &Vec<f64>, n: usize, h: f64, epsilon: f64) -> Vec<f64>
{
    let dx = f_x(x,y,n,epsilon);
    let answer: Vec<f64> = (0..n).map(|i| x[i] + h * dx[i]).collect();
    answer
}

fn step_y(x: &Vec<f64>, y: &Vec<f64>, n: usize, h: f64, epsilon: f64) -> Vec<f64>
{
    let dy = f_y(x,y,n,epsilon);
    let answer: Vec<f64> = (0..n).map(|i| y[i] + h * dy[i]).collect();
    answer
}

fn rk4_x(x: &[f64], y: &[f64], n: usize, h: f64, epsilon: f64) -> Vec<f64>
{
    let k1 : Vec<f64> = f_x(x,y,n,epsilon).iter().map(|u| u * h).collect();
    let k0_k12: Vec<f64> = x.iter().zip(k1.iter()).map(|(a,b)| a + b / 2.0).collect();
    let k2 : Vec<f64> = f_x(&k0_k12,y,n,epsilon).iter().map(|u| u * h).collect();
    let k0_k22 : Vec<f64> = x.iter().zip(k2.iter()).map(|(a,b)| a + b / 2.0).collect();
    let k3 : Vec<f64> = f_x(&k0_k22,y,n,epsilon).iter().map(|u| u * h).collect();
    let k0_k3 : Vec<f64> = x.iter().zip(k3.iter()).map(|(a,b)| a + b).collect();
    let k4 : Vec<f64> = f_x(&k0_k3, y,n,epsilon).iter().map(|u| u * h).collect();

    k1.iter().zip(k2.iter()).map(|(a,b)| a + 2.0 * b)
        .zip(k3.iter()).map(|(a,b)| a + 2.0 * b)
        .zip(k4.iter()).map(|(a,b)| (a + b) / 6.0) 
        .zip(x.iter()).map(|(a,b)| a + b).collect()
}

fn rk4_y(x: &[f64], y: &[f64], n: usize, h: f64, epsilon: f64) -> Vec<f64>
{
    let k1 : Vec<f64> = f_y(x,y,n,epsilon).iter().map(|u| u * h).collect();
    let k0_k12: Vec<f64> = y.iter().zip(k1.iter()).map(|(a,b)| a + b / 2.0).collect();
    let k2 : Vec<f64> = f_y(x,&k0_k12,n,epsilon).iter().map(|u| u * h).collect();
    let k0_k22 : Vec<f64> = y.iter().zip(k2.iter()).map(|(a,b)| a + b / 2.0).collect();
    let k3 : Vec<f64> = f_y(x,&k0_k22,n,epsilon).iter().map(|u| u * h).collect();
    let k0_k3 : Vec<f64> = y.iter().zip(k3.iter()).map(|(a,b)| a + b).collect();
    let k4 : Vec<f64> = f_y(x,&k0_k3,n,epsilon).iter().map(|u| u * h).collect();

    k1.iter().zip(k2.iter()).map(|(a,b)| a + 2.0 * b)
        .zip(k3.iter()).map(|(a,b)| a + 2.0 * b)
        .zip(k4.iter()).map(|(a,b)| (a + b) / 6.0) 
        .zip(y.iter()).map(|(a,b)| a + b).collect()
}

fn init_x(n: usize) -> Vec<f64> {
    let d_xi = 1.0 / (0.0 + (n as f64));
    let mut xi = 0.0;

    let mut answer : Vec<f64> = Vec::new();

    for _i in 1..=n 
    {
        xi = xi + d_xi;
        answer.push(xi + 0.01 * (2.0 * 3.0 * std::f64::consts::PI * xi).sin() );
    }
    answer
}

fn init_y(n: usize) -> Vec<f64> {
    let d_xi = 1.0 / (0.0 + (n as f64));
    let mut xi = 0.0;

    let mut answer : Vec<f64> = Vec::new();

    for _i in 1..=n 
    {
        xi = xi + d_xi;
        answer.push(-0.01 * (2.0 * 3.0 * std::f64::consts::PI * xi).sin() );
    }
    answer
}

fn f_x(x: &[f64], y: &[f64], n: usize, epsilon: f64) -> Vec<f64> 
{
    let mut answer: Vec<f64> = vec![];
    answer.resize(n, 0.0);

    for i in 0..n
    {
        let mut accum: f64 = 0.0;
        for j in 0..n 
        {
            if i == j {
                continue;
            }
            let mut temp: f64; 
            temp = ( (y[i] - y[j]) * 2.0 * std::f64::consts::PI).sinh();
            temp /= 
                ((y[i] - y[j]) * 2.0 * std::f64::consts::PI).cosh() 
                - ((x[i] - x[j]) * 2.0 * std::f64::consts::PI).cos() 
                + epsilon * epsilon
            ;
            accum += temp;
        }
        accum /= 2.0 * (n as f64);
        accum *= -1.0;
        answer[i] = accum;
    }
    answer
}

fn f_y(x: &[f64], y: &[f64], n: usize, epsilon: f64) -> Vec<f64> 
{
    let mut answer: Vec<f64> = vec![];
    answer.resize(n, 0.0);

    for i in 0..n
    {
        let mut accum: f64 = 0.0;
        for j in 0..n 
        {
            if i == j {
                continue;
            }
            let mut temp: f64; 
            temp = ( (x[i] - x[j]) * 2.0 * std::f64::consts::PI).sin();
            temp /= 
                ((y[i] - y[j]) * 2.0 * std::f64::consts::PI).cosh() 
                - ((x[i] - x[j]) * 2.0 * std::f64::consts::PI).cos() 
                + epsilon * epsilon
            ;
            accum += temp;
        }
        accum /= 2.0 * (n as f64);
        accum *= 1.0;
        answer[i] = accum;
    }
    answer
}