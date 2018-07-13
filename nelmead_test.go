package nelmead

import (
	"testing"
)

func Rosenbrock(x []float64) float64 {
	const a = 1.0
	const b = 100.0
	return (a-x[0])*(a-x[0]) + b*(x[1]-x[0]*x[0])*(x[1]-x[0]*x[0])
}

func TestSimplexValues_Centroid(t *testing.T) {
	s := MakeSimplexValues([]float64{0.0, 0.0, 0.0}, 3.0)
	x := s.Centroid()
	if !(x[0] == 1. && x[1] == 1.0 && x[2] == 0.0) {
		t.Errorf("expected [1. 1. 0.] got: %v", x)
	}
}

func TestNelderMeadOptimizer_Optimize(t *testing.T) {
	const tol = 1.E-9
	start := []float64{-1.0, 1.0}
	step := 0.5
	nm := MakeNelderMeadOptimizer(Rosenbrock, 2)
	xval, opt, converged := nm.Optimize(200, start, step)
	if xval[0]-1.0 > tol || xval[1]-1.0 > tol || !converged {
		t.Errorf("%v %v %v", xval, opt, converged)
	}
}
