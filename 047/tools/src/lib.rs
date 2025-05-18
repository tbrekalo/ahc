#![allow(non_snake_case)]

use itertools::Itertools;
use proconio::{input, marker::Chars};
use rand::prelude::*;
use std::{collections::HashMap, ops::RangeBounds};
use svg::node::element::{Circle, Definitions, Group, Marker, Path, Rectangle, Style, Title};

pub trait SetMinMax {
    fn setmin(&mut self, v: Self) -> bool;
    fn setmax(&mut self, v: Self) -> bool;
}
impl<T> SetMinMax for T
where
    T: PartialOrd,
{
    fn setmin(&mut self, v: T) -> bool {
        *self > v && {
            *self = v;
            true
        }
    }
    fn setmax(&mut self, v: T) -> bool {
        *self < v && {
            *self = v;
            true
        }
    }
}

#[macro_export]
macro_rules! mat {
    ($($e:expr),*) => { Vec::from(vec![$($e),*]) };
    ($($e:expr,)*) => { Vec::from(vec![$($e),*]) };
    ($e:expr; $d:expr) => { Vec::from(vec![$e; $d]) };
    ($e:expr; $d:expr $(; $ds:expr)+) => { Vec::from(vec![mat![$e $(; $ds)*]; $d]) };
}

#[derive(Clone, Debug)]
pub struct Input {
    pub N: usize,
    pub M: usize,
    pub L: usize,
    pub S: Vec<Vec<char>>,
    pub P: Vec<i64>,
}

impl std::fmt::Display for Input {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "{} {} {}", self.N, self.M, self.L)?;
        for i in 0..self.N {
            writeln!(f, "{} {}", self.S[i].iter().collect::<String>(), self.P[i])?;
        }
        Ok(())
    }
}

#[derive(Clone, Debug)]
pub struct Output {
    pub out: Vec<SingleOutput>,
}

#[derive(Clone, Debug)]
pub struct SingleOutput {
    pub C: Vec<char>,
    pub A: Vec<Vec<i32>>,
}

pub fn parse_input(f: &str) -> Input {
    let f = proconio::source::once::OnceSource::from(f);
    input! {
        from f,
        N: usize, M: usize, L: usize,
        words_and_scores: [(Chars, i64); N]
    }

    let (S, P) = words_and_scores.into_iter().unzip();

    Input { N, M, L, S, P }
}

pub fn read<T: Copy + PartialOrd + std::fmt::Display + std::str::FromStr, R: RangeBounds<T>>(
    token: Option<&str>,
    range: R,
) -> Result<T, String> {
    if let Some(v) = token {
        if let Ok(v) = v.parse::<T>() {
            if !range.contains(&v) {
                Err(format!("Out of range: {}", v))
            } else {
                Ok(v)
            }
        } else {
            Err(format!("Parse error: {}", v))
        }
    } else {
        Err("Unexpected EOF".to_owned())
    }
}

pub fn parse_output(input: &Input, f: &str) -> Result<Output, String> {
    let mut f = f.split_whitespace().peekable();
    let mut out = vec![];
    while f.peek().is_some() {
        let mut C = vec!['.'; input.M];
        let mut A = mat![0; input.M; input.M];
        for i in 0..input.M {
            C[i] = read(f.next(), 'a'..='f')?;
            for j in 0..input.M {
                A[i][j] = read(f.next(), 0..=100)?;
            }
        }
        out.push(SingleOutput { C, A });
    }
    if out.len() == 0 {
        return Err("Empty Output".to_owned());
    }
    Ok(Output { out })
}

pub fn compute_score(input: &Input, out: &Output) -> (i64, String) {
    let (mut score, err, _) = compute_score_details(input, &out.out.last().unwrap());
    if err.len() > 0 {
        score = 0;
    }
    (score, err)
}

pub fn compute_score_details(input: &Input, out: &SingleOutput) -> (i64, String, Vec<f64>) {
    for i in 0..input.M {
        let sum = out.A[i].iter().sum::<i32>();
        if sum != 100 {
            return (
                0,
                format!("Sum of A[{}] is not 100 ({})", i, sum),
                vec![0.0; input.N],
            );
        }
    }
    let mut total_score = 0.0;
    let mut probs = vec![];
    for (word, &score) in input.S.iter().zip(&input.P) {
        let prob = compute_word_probability(word, input.L, &out.C, &out.A);
        probs.push(prob);
        total_score += prob * score as f64;
    }
    (total_score.round() as i64, String::new(), probs)
}

fn compute_word_probability(word: &Vec<char>, L: usize, C: &Vec<char>, A: &Vec<Vec<i32>>) -> f64 {
    let M = C.len();
    // Enumerate the automaton states (matching length, vertex)
    let mut n = 0;
    let mut states = HashMap::new();
    for j in 0..M {
        states.insert((0, j), n);
        n += 1;
        for i in 0..word.len() - 1 {
            if word[i] == C[j] {
                states.insert((i + 1, j), n);
                n += 1;
            }
        }
    }
    // Construct the state transition matrix
    let mut X = vec![0.0; n * n];
    for (&(len, u), &j) in &states {
        for v in 0..M {
            let mut next = word[..len].to_vec();
            next.push(C[v]);
            let mut s = 0;
            while next[s..] != word[..next.len() - s] {
                s += 1;
            }
            if next.len() - s != word.len() {
                let i = states[&(next.len() - s, v)];
                X[i * n + j] += A[u][v] as f64 / 100.0;
            }
        }
    }
    // Compute Y=X^(L-1)
    let mut power = L - 1;
    let mut Y = vec![0.0; n * n];
    for i in 0..n {
        Y[i * n + i] = 1.0;
    }
    while power > 0 {
        if power & 1 != 0 {
            Y = mul(&Y, &X, n);
        }
        X = mul(&X, &X, n);
        power >>= 1;
    }
    // Compute the probability
    let init = if C[0] == word[0] {
        states[&(1, 0)]
    } else {
        states[&(0, 0)]
    };
    let mut ret = 1.0;
    for i in 0..n {
        ret -= Y[i * n + init];
    }
    ret.clamp(0.0, 1.0)
}

fn mul(a: &Vec<f64>, b: &Vec<f64>, n: usize) -> Vec<f64> {
    let mut c = vec![0.0; n * n];
    for i in 0..n {
        for k in 0..n {
            for j in 0..n {
                c[i * n + j] += a[i * n + k] * b[k * n + j];
            }
        }
    }
    c
}

pub fn gen(seed: u64) -> Input {
    let mut rng = rand_chacha::ChaCha20Rng::seed_from_u64(seed);
    let N = 36;
    let M = 12;
    let L = 1_000_000;
    let mut S = vec![];
    let mut P = vec![];
    let chars = vec!['a', 'b', 'c', 'd', 'e', 'f'];
    for _ in 0..N {
        loop {
            let len = rng.gen_range(6..=12);
            let word = (0..len).map(|_| *chars.choose(&mut rng).unwrap()).collect();
            if S.contains(&word) {
                continue;
            }
            let max_score = (1.5f64.powi(2 * len as i32)) as i64;
            let score = rng.gen_range(1..=max_score);
            S.push(word);
            P.push(score);
            break;
        }
    }
    Input { N, M, L, S, P }
}

/// 0 <= val <= 1
pub fn color(mut val: f64) -> String {
    val.setmin(1.0);
    val.setmax(0.0);
    let (r, g, b) = if val < 0.5 {
        let x = val * 2.0;
        (
            30. * (1.0 - x) + 144. * x,
            144. * (1.0 - x) + 255. * x,
            255. * (1.0 - x) + 30. * x,
        )
    } else {
        let x = val * 2.0 - 1.0;
        (
            144. * (1.0 - x) + 255. * x,
            255. * (1.0 - x) + 30. * x,
            30. * (1.0 - x) + 70. * x,
        )
    };
    format!(
        "#{:02x}{:02x}{:02x}c0",
        r.round() as i32,
        g.round() as i32,
        b.round() as i32
    )
}

pub fn rect(x: usize, y: usize, w: usize, h: usize, fill: &str) -> Rectangle {
    Rectangle::new()
        .set("x", x)
        .set("y", y)
        .set("width", w)
        .set("height", h)
        .set("fill", fill)
}

pub fn group(title: String) -> Group {
    Group::new().add(Title::new(title))
}

pub fn vis_default(input: &Input, out: &Output) -> (i64, String, String) {
    let (mut score, err, svg) = vis(input, &out.out.last().unwrap());
    if err.len() > 0 {
        score = 0;
    }
    (score, err, svg)
}

pub fn vis(input: &Input, out: &SingleOutput) -> (i64, String, String) {
    let W = 1000;
    let H = 500;
    let (score, err, probs) = compute_score_details(input, out);
    let mut doc = svg::Document::new()
        .set("id", "vis")
        .set("viewBox", (-5, -5, W + 10, H + 10))
        .set("width", W + 10)
        .set("height", H + 10)
        .set("style", "background-color:white");
    doc = doc.add(Style::new("text { font-family:  monospace; }"));
    let arrow = Marker::new()
        .set("id", "arrow")
        .set("markerUnits", "userSpaceOnUse")
        .set("viewBox", "0 0 8 8")
        .set("refX", 8)
        .set("refY", 4)
        .set("markerWidth", 8)
        .set("markerHeight", 8)
        .set("orient", "auto")
        .add(
            Path::new()
                .set("d", "M 0 0 L 8 4 L 0 8 z")
                .set("fill", "context-stroke"),
        );
    doc = doc.add(Definitions::new().add(arrow));
    let r = H as f64 * 0.4;
    let c = H as f64 * 0.5;
    let cs = (0..input.M)
        .map(|i| {
            let rad = std::f64::consts::PI * (2.0 * i as f64 / input.M as f64 - 0.5);
            let x = c + r * rad.cos();
            let y = c + r * rad.sin();
            (x, y)
        })
        .collect_vec();
    for i in 0..input.M {
        for j in 0..input.M {
            if out.A[i][j] > 0 {
                if i != j {
                    let (x1, y1) = cs[i];
                    let (x2, y2) = cs[j];

                    let dx = x2 - x1;
                    let dy = y2 - y1;
                    let len = (dx * dx + dy * dy).sqrt();
                    let (px, py) = (-dy / len, dx / len);
                    let offset_amount = if out.A[j][i] > 0 { 20.0 } else { 0.0 };
                    let mx = (x1 + x2) * 0.5;
                    let my = (y1 + y2) * 0.5;
                    let cx = mx + px * offset_amount;
                    let cy = my + py * offset_amount;

                    let r = 8.0;
                    let d0 = ((cx - x1).hypot(cy - y1)).max(1e-6);
                    let d2 = ((x2 - cx).hypot(y2 - cy)).max(1e-6);
                    let t0 = (r / d0).clamp(0.0, 1.0);
                    let t1 = (1.0 - r / d2).clamp(0.0, 1.0);

                    let a0x = x1 + (cx - x1) * t0;
                    let a0y = y1 + (cy - y1) * t0;
                    let b0x = cx + (x2 - cx) * t0;
                    let b0y = cy + (y2 - cy) * t0;
                    let c0x = a0x + (b0x - a0x) * t0;
                    let c0y = a0y + (b0y - a0y) * t0;

                    let a1x = x1 + (cx - x1) * t1;
                    let a1y = y1 + (cy - y1) * t1;
                    let b1x = cx + (x2 - cx) * t1;
                    let b1y = cy + (y2 - cy) * t1;
                    let c1x = a1x + (b1x - a1x) * t1;
                    let c1y = a1y + (b1y - a1y) * t1;

                    let ratio = if (1.0 - t0) != 0.0 {
                        (t1 - t0) / (1.0 - t0)
                    } else {
                        0.0
                    };
                    let d0x = a0x + (b0x - a0x) * ratio;
                    let d0y = a0y + (b0y - a0y) * ratio;

                    doc = doc.add(
                        group(format!("A[{}][{}] = {}", i, j, out.A[i][j])).add(
                            Path::new()
                                .set(
                                    "d",
                                    format!("M {} {} Q {} {} {} {}", c0x, c0y, d0x, d0y, c1x, c1y),
                                )
                                .set("marker-end", "url(#arrow)")
                                // .set("stroke", color(out.A[i][j] as f64 / 100.0))
                                .set(
                                    "stroke",
                                    format!(
                                        "#000000{:02x}",
                                        (30 + 255 * out.A[i][j] / 100).min(255)
                                    ),
                                )
                                .set("stroke-width", 1.0 + out.A[i][j] as f64 / 40.0)
                                .set("fill", "none"),
                        ),
                    );
                } else {
                    let (x, y) = cs[i];
                    let r0 = 15.0;
                    let loop_r = 20.0;
                    let (cx0, cy0) = (H as f64 * 0.5, H as f64 * 0.5);
                    let dx = x - cx0;
                    let dy = y - cy0;
                    let d = (dx * dx + dy * dy).sqrt().max(1e-6);
                    let (ux, uy) = (dx / d, dy / d);
                    let sx = x + ux * r0;
                    let sy = y + uy * r0;
                    let (px, py) = (-uy, ux);
                    let c1x = sx + ux * loop_r + px * (loop_r);
                    let c1y = sy + uy * loop_r + py * (loop_r);
                    let c2x = sx + ux * loop_r - px * (loop_r);
                    let c2y = sy + uy * loop_r - py * (loop_r);
                    doc = doc.add(
                        group(format!("A[{}][{}] = {}", i, i, out.A[i][i])).add(
                            Path::new()
                                .set(
                                    "d",
                                    format!(
                                        "M {} {} C {} {} {} {} {} {}",
                                        sx, sy, c1x, c1y, c2x, c2y, sx, sy
                                    ),
                                )
                                .set("marker-end", "url(#arrow)")
                                // .set("stroke", color(out.A[i][i] as f64 / 100.0))
                                .set(
                                    "stroke",
                                    format!(
                                        "#000000{:02x}",
                                        (30 + 255 * out.A[i][j] / 100).min(255)
                                    ),
                                )
                                .set("stroke-width", 1.0 + out.A[i][j] as f64 / 40.0)
                                .set("fill", "none"),
                        ),
                    );
                }
            }
        }
    }

    for i in 0..input.M {
        let (x, y) = cs[i];
        doc = doc.add(
            group(format!("State {}", i))
                .add(
                    Circle::new()
                        .set("cx", x)
                        .set("cy", y)
                        .set("r", 15)
                        .set("fill", "white")
                        .set("stroke", "black")
                        .set("stroke-width", 1),
                )
                .add(
                    svg::node::element::Text::new(out.C[i])
                        .set("x", x)
                        .set("y", y)
                        .set("font-size", 20)
                        .set("fill", "black")
                        .set("text-anchor", "middle")
                        .set("dominant-baseline", "middle"),
                ),
        );
    }

    let D = H / (input.N / 2);
    for i in 0..input.N {
        let d = i / (input.N / 2);
        doc = doc.add(
            group(format!("prob[{}] = {:.10}", i, probs[i]))
                .add(rect(
                    H + (W - H) / 2 * d,
                    D * (i % (input.N / 2)),
                    (W - H) / 2 - 5,
                    D,
                    "white",
                ))
                .add(rect(
                    H + (W - H) / 2 * d,
                    D * (i % (input.N / 2)),
                    (((W - H) / 2 - 5) as f64 * probs[i]).round() as usize,
                    D,
                    &color(probs[i]),
                ))
                .add(
                    svg::node::element::Text::new(format!("{}", i))
                        .set("x", H + (W - H) / 2 * d + 5)
                        .set("y", D * (i % (input.N / 2)) + D / 2)
                        .set("font-size", 20)
                        .set("fill", "black")
                        .set("text-anchor", "start")
                        .set("dominant-baseline", "central"),
                )
                .add(
                    svg::node::element::Text::new(format!("{}", input.P[i]))
                        .set("x", H + (W - H) / 2 * d + 85)
                        .set("y", D * (i % (input.N / 2)) + D / 2)
                        .set("font-size", 20)
                        .set("fill", "black")
                        .set("text-anchor", "end")
                        .set("dominant-baseline", "central"),
                )
                .add(
                    svg::node::element::Text::new(format!(
                        "{}",
                        input.S[i].iter().collect::<String>()
                    ))
                    .set("x", H + (W - H) / 2 * d + 100)
                    .set("y", D * (i % (input.N / 2)) + D / 2)
                    .set("font-size", 20)
                    .set("fill", "black")
                    .set("text-anchor", "start")
                    .set("dominant-baseline", "central"),
                ),
        );
    }
    for i in 0..input.N {
        let d = i / (input.N / 2);
        doc = doc.add(
            rect(
                H + (W - H) / 2 * d,
                D * (i % (input.N / 2)),
                (W - H) / 2 - 5,
                D,
                "none",
            )
            .set("stroke", "black")
            .set("stroke-width", 1),
        );
    }

    (score, err, doc.to_string())
}
