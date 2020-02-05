#[macro_use]
mod kattis;

use std::collections::BinaryHeap;
use std::i32;

struct FlowNetwork {
    v: usize,
    s: usize,
    t: usize,
    c: Vec<i32>,
    f: Vec<i32>,
    adj: Vec<Vec<usize>>,
    x: Vec<i32>,
    l: Vec<usize>,
    pq: BinaryHeap<usize>,
}

impl FlowNetwork {
    fn new(v: usize, s: usize, t: usize) -> Self {
        Self {
            v,
            s,
            t,
            c: vec![0; v * v],
            f: vec![0; v * v],
            adj: vec![Vec::new(); v],
            x: vec![0; v],
            l: vec![0; v],
            pq: BinaryHeap::new(),
        }
    }

    fn insert_edge(&mut self, a: usize, b: usize, c: u32) {
        self.set_c(a, b, c);
        self.adj[a].push(b);
        self.adj[b].push(a);
    }

    fn init(&mut self) {
        self.l[self.s] = self.v;
        for i in 0..self.adj[self.s].len() {
            let b = self.adj[self.s][i];
            self.mod_f(self.s, b, self.c(self.s, b));
        }
    }

    fn calculate_flow(&mut self) {
        while let Some(a) = self.pq.pop() {
            self.discharge(a);
        }
    }

    fn total_flow(&self) -> usize {
        let mut flow = 0;
        for i in 0..self.adj[self.t].len() {
            let a = self.adj[self.t][i];
            if self.f(a, self.t) > 0 {
                flow += self.f(a, self.t) as usize;
            }
        }
        flow
    }

    fn discharge(&mut self, a: usize) {
        let mut i = 0;
        while self.x[a] > 0 {
            if i == self.adj[a].len() {
                self.relabel(a);
                i = 0;
            } else {
                let b = self.adj[a][i];
                if self.l[a] == self.l[b] + 1 {
                    self.push(a, b);
                }
                i += 1;
            }
        }
    }

    fn push(&mut self, a: usize, b: usize) {
        let delta = i32::min(self.x[a], self.cf(a, b));
        self.mod_f(a, b, delta)
    }

    fn relabel(&mut self, a: usize) {
        let mut min = 2 * self.v;
        for i in 0..self.adj[a].len() {
            let b = self.adj[a][i];
            if self.cf(a, b) > 0 && self.l[b] < min {
                min = self.l[b];
            }
        }
        self.l[a] = min + 1;
    }

    fn set_c(&mut self, a: usize, b: usize, c: u32) {
        self.c[b * self.v + a] = c as i32;
    }

    fn mod_f(&mut self, a: usize, b: usize, delta: i32) {
        self.f[b * self.v + a] += delta;
        self.f[a * self.v + b] -= delta;
        self.x[a] -= delta;
        self.x[b] += delta;
        if b != self.t && self.x[b] <= delta {
            self.pq.push(b);
        }
    }

    fn c(&self, a: usize, b: usize) -> i32 {
        self.c[b * self.v + a]
    }

    fn f(&self, a: usize, b: usize) -> i32 {
        self.f[b * self.v + a]
    }

    fn cf(&self, a: usize, b: usize) -> i32 {
        self.c[b * self.v + a] - self.f[b * self.v + a]
    }
}

#[allow(clippy::many_single_char_names)]
fn main() {
    let (v, s, t, e): (usize, usize, usize, usize);
    scanln!(v);
    scanln!(s, t);
    scanln!(e);
    let mut f_net = FlowNetwork::new(v, s - 1, t - 1);
    for _ in 0..e {
        let (a, b, c): (usize, usize, u32);
        scanln!(a, b, c);
        if b != s && a != t {
            f_net.insert_edge(a - 1, b - 1, c);
        }
    }
    f_net.init();
    f_net.calculate_flow();
    println!("{}", f_net.v);
    println!("{} {} {}", f_net.s + 1, f_net.t + 1, f_net.total_flow());
    let mut edges = Vec::new();
    for a in 0..v {
        for b in 0..v {
            if f_net.f(a, b) > 0 {
                edges.push((a, b, f_net.f(a, b)));
            }
        }
    }
    println!("{}", edges.len());
    for (a, b, f) in edges {
        println!("{} {} {}", a + 1, b + 1, f);
    }
}
