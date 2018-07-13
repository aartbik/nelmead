# Nelder-Mead Optimization in Go 

Package `nelmead` is s a small implementation of the Nelder-Mead optimization
algorithm in Go with minimal dependencies. It strictly follows the algorithm as
described in  https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method. 

This is also my first contribution to GitHub, and a bit of a trial for me how
this works.

It requires a yield or optimizing function, and the number of parameters to be
optimized (dim)

```go
	func MakeNelderMeadOptimizer(f func([]float64) float64, dim int) *NelderMeadOptimizer 
```

Default values for the optimizer can be used, but if you want to tweak the
performance you can adjust various parameter in the `struct
NelderMeadOptimizer`.

```go
type NelderMeadOptimizer struct {

    // ConvTreshold is the minimimal difference between two successive trials
    // for an optimization to be classified as an improvement.
    ConvTreshold float64

    // ConvLimit is the number of iterations with no significant improvement in
    // the optimizing function for the solution to be classified as converged.
    ConvLimit int

	// These are the parameters used for the reflection, expansion,
	// contraction, and shrink stages in the optimzation process
    Alpha     float64
    Gamma     float64
    Rho       float64
    Sigma     float64
}
```

The actual optimization trial is started by using the Optimize method
associated with an instance of the `NelderMeadOptimizer`:
```go
func (nmo *NelderMeadOptimizer) Optimize(max_iter int, start []float64, step float64) ([]float64, float64, bool) 
```
As input it requires the maximum number of iterations `max_iter`, and
a `start`-value for the parameters of the function, and a `step`-size. 

It returns  the solution vector, the value of the minimum found for the
solution, and a boolean indicating if the solution has converged or not.

The solution is strongly dependent on the `start` and `step` values and might
require tweaking to get the best solution. For more detailed information see
godoc. Also see the `nelmead_test.go` file as example how to use this function.

# License

Copyright 2018 Gerard Harbers

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this file except in compliance with the License. You may obtain a copy of the
License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.


