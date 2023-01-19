---
title: "RK4 in different languages"
draft: true
---

||Julia | C ++ | Python | 
|----|------|------|--------|
|**speeds**|1.2|3.4|5.6|

## `Julia`
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

    plot(output.times, output.solutions, label="RK4")
    plot!(output.times, map(exact_func, output.times), ls=:dash, lw=2, label="exact")
    title!("RK4 in julia")
    xlabel!("time")
    ylabel!(L"df/dt")
    savefig("rk4_julia.png")
end
main()

```
![julia plot](/rk4_julia.png#center)
## `C++`
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

void RungeKutta::rk4Solver(float (*func)(float, float), float initialSolution) {
  vector<float> times(this->numSteps + 1);
  times = this->createTimes(0.0);
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
  RungeKutta rkObj = RungeKutta(100, 0.1);

  auto t1 = chrono::high_resolution_clock::now();
  rkObj.rk4Solver(&exampleFunction, 1.0);
  auto t2 = chrono::high_resolution_clock::now();
  cout << chrono::duration_cast<chrono::nanoseconds>(t2-t1).count() * pow(10, -9);
  return 0;
} 

```
![cpp plot](/rk4_cpp.jpg#center)
## Python
```
import numpy as np
import matplotlib.pyplot as plt

class Solve():
    def __init__(self, num_steps: int, step_size: float, constant: float):
        self.constant = constant
        self.num_steps = num_steps
        self.step_size = step_size

    def create_times(self, initial_time):
        T = np.array([initial_time + n * self.step_size for n in range(self.num_steps + 1)])
        return T

    def runge_kutta(self, f, times: np.array, s_0=1):

        S = np.zeros(self.num_steps + 1)
            
        S[0] = s_0

        for n in range(self.num_steps):
            k1 = self.step_size * f(times[n], S[n])
            k2 = self.step_size * f(times[n] + self.step_size / 2, S[n] + k1 / 2)
            k3 = self.step_size * f(times[n] + self.step_size / 2, S[n] + k2 / 2)
            k4 = self.step_size * f(times[n] + self.step_size, S[n] + k3)
            S[n + 1] = S[n] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
                
        return S

    def test_func(self, t, x):
        y = - self.constant * x
        return y

    def exact_solution(self, t):
        y = - np.exp(- self.constant * t)
        return y

s = Solve(1000, 0.01, 10)
times = s.create_times(0)
sols = s.runge_kutta(s.test_func, times)
plt.plot(times, sols)
plt.show()
```