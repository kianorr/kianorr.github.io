---
title: "RK4 speeds in python, julia, and c++"
draft: false
date: 2023-01-20
lastmod: 2023-02-12
math: true
tags: ["python", "julia", "c++"]
summary: "I wanted look more into Julia because it sounded like python but faster and
better. So, I compared Julia to Python, and also threw C++ in there 
(do I regret that? Maybe)"
tocOpen: true
---

| Julia | C ++ | Python | 
| :------: | :------: | :--------: |
| 22.555 ± 34.788 μs | 14.162 μs | 2420 ± 197 μs |

## intro
Fourth order Runge Kutta (RK4) is a computational method to solve DEs
I wanted look more into Julia because it sounded like python but faster and
better. So, I tested speeds between Julia and Python, and also threw C++ in there 
(do I regret that? Maybe), by using fourth order Runge Kutta (RK4). RK4 is a computational method
to solve differential equations, which in my case was the simple

\begin{align}
\frac{\text{d}x}{\text{d}t} = - 10x.
\end{align}

Julia is actually quite close to C++, which isn't too surprising, but it's nice to see.
Now everyone just needs to actually use Julia so there's support behind it :p

## `Julia`
I tried learning how to make my code fast as best I could, but I know there's
still a lot of room for improvement. I'm still trying to figure out the best
way to access the variables in structs. I've just been making them input parameters
in functions.

```
using Plots, LaTeXStrings, BenchmarkTools
theme(:dark)

struct InitialConditions
    initial_solution::Float64
    initial_time::Float64
    numberof_steps::Int
    step_size::Float64
end

mutable struct Output{T}
    times::Vector{T}
    solutions::Vector{T}
end

function make_times(ic::InitialConditions, out::Output)
    out.times::Vector{Float64} = map(n -> ic.initial_time + n * ic.step_size,
    				     range(1, length=ic.numberof_steps + 1))
end

function rk4(f::Function, ic::InitialConditions, out::Output)
    out.solutions[1] = ic.initial_solution

    for n in range(1, length=ic.numberof_steps)
        k1::Float64 = ic.step_size .* f(out.times[n], out.solutions[n])
        k2::Float64 = ic.step_size .* f(out.times[n] .+ ic.step_size / 2, out.solutions[n] .+ k1 ./ 2)
        k3::Float64 = ic.step_size .* f(out.times[n] .+ ic.step_size / 2, out.solutions[n] .+ k2 ./ 2)
        k4::Float64 = ic.step_size .* f(out.times[n] .+ ic.step_size, out.solutions[n] .+ k3)
        out.solutions[n + 1] = out.solutions[n] .+ (k1 .+ 2 * k2 .+ 2 .* k3 .+ k4) ./ 6
    end
end

test_func(t, x) = - 10 * x
exact_func(t) = exp(- 10 * t)
function main()
    input::InitialConditions = InitialConditions(1.0, 0.0, 1000, 0.001)
    output::Output = Output(Vector{Float64}(undef, input.numberof_steps + 1), 
                            Vector{Float64}(undef, input.numberof_steps + 1))

    make_times(input, output)
    rk4(test_func, input, output)
end

function plot_rk4()
    output = main()
    plot(output.times, output.solutions, label="RK4")
    plot!(output.times, map(exact_func, output.times), ls=:dash, lw=2, label="exact")
    title!("RK4 in julia")
    xlabel!("time")
    ylabel!(L"df/dt")
    savefig("rk4_julia.png")
end

@benchmark main()
plot_rk4()

```
Julia has this really pretty benchmarking macro that will print this in the
Julia REPL:
![julia benchmark](/julia_benchmark.png#center)

and the output for the plot is
![julia plot](/rk4_julia.png#center)

## `C++`
oh boy

```
#include <iostream>
#include <chrono>
#include <math.h>
#include <vector>
using namespace std;

class RungeKutta {
  public:
    void rk4Solver(float (*f)(float, float), float s);

    int numSteps;
    float stepSize;

    RungeKutta(int n, float h) {
      numSteps = n;
      stepSize = h;
    }
  
  private:
    vector<float> createTimes(float t);
};

vector<float> RungeKutta::createTimes(float initialTime) {
  vector<float> times(this->numSteps + 1);
  for (int i; i <= this->numSteps + 1; i ++) {
    times[i] = i;
  }
  return times;
}

void RungeKutta::rk4Solver(float (*func)(float, float), float initialSolution, float initialTime) {
  vector<float> times(this->numSteps + 1);
  times = this->createTimes(initialTime);
  float solutions[this->numSteps + 1];

  float k1;
  float k2;
  float k3;
  float k4;

  solutions[0] = initialSolution;
  for (int n; n <= this->numSteps; n ++) {
    k1 = this->stepSize * func(times[n], solutions[n]);
    k2 = this->stepSize * func(times[n] + this->stepSize / 2, solutions[n] + k1 / 2);
    k3 = this->stepSize * func(times[n] + this->stepSize / 2, solutions[n] + k2 / 2);
    k4 = this->stepSize * func(times[n] + this->stepSize, solutions[n] + k3);
    solutions[n + 1] = solutions[n] + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
  }
}

float exampleFunction(float t, float x) {
  float y = - 10 * x;
  return y;
}

float exactSolution(float t) {
  float y = - exp(- 10 * t);
  return y;
}

int main() {
  RungeKutta rkObj = RungeKutta(1000, 0.001);

  auto t1 = chrono::high_resolution_clock::now();
  rkObj.rk4Solver(&exampleFunction, 1.0, 0.0);
  auto t2 = chrono::high_resolution_clock::now();
  cout << chrono::duration_cast<chrono::nanoseconds>(t2-t1).count() * pow(10, -9);
  return 0;
} 

```
I had to use `ROOT` to get plots lol. I like the aesthetic of `ROOT` plots
though.

![cpp plot](/rk4_cpp.jpg#center)

## `python`
ahh python. This is one of the first times I'm using `@staticmethod` so I think
that helped with the speed a bit.
```
import numpy as np
import matplotlib.pyplot as plt

class Solve():
    def __init__(self, initial_solution: float, initial_time: float,
                 num_steps: int, step_size: float):

        self.initial_time = initial_time
        self.initial_solution = initial_solution

        self.numberof_steps = num_steps
        self.step_size = step_size

    def create_times(self):
        times = np.array([self.initial_time + n * self.step_size for n in range(self.numberof_steps + 1)])
        return times

    def find_solutions_via_rk4(self, f, times: np.array):

        solutions = np.zeros(self.numberof_steps + 1)
            
        solutions[0] = self.initial_solution

        for n in range(self.numberof_steps):
            k1 = self.step_size * f(times[n], solutions[n])
            k2 = self.step_size * f(times[n] + self.step_size / 2, solutions[n] + k1 / 2)
            k3 = self.step_size * f(times[n] + self.step_size / 2, solutions[n] + k2 / 2)
            k4 = self.step_size * f(times[n] + self.step_size, solutions[n] + k3)
            solutions[n + 1] = solutions[n] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
                
        return solutions

    @staticmethod
    def test_func(t, x):
        y = - 10 * x
        return y

    @staticmethod
    def exact_solution(t):
        y = - np.exp(- 10 * t)
        return y
```

and then for speed tests, I used `%%timeit` in jupyter notebook:
```
from rk4 import Solve
import matplotlib.pyplot as plt
plt.style.use("ggplot")
```
```
%%timeit
s = Solve(1.0, 0.0, 1000, 0.001)
times = s.create_times()
sols = s.find_solutions_via_rk4(s.test_func, times)
```
Output:
`2.42 ms ± 197 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)`
```
fig, ax = plt.subplots(figsize=(8, 6))
plt.plot(times, sols)
plt.ylabel("$df / dt$", fontsize=15)
plt.xlabel("time", fontsize=15)
plt.title("RK4 in Python", fontsize=17)
plt.tight_layout()
plt.savefig("rk4_python.png")
plt.show()
```
![rk4 in python](/rk4_python.png#center)