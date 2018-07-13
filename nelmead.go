// Copyright 2018 Gerard Harbers
//
// Licensed under the Apache License, Version 2.0 (the "License"); you may
// not use this file except in compliance with the License. You may obtain a
//
// copy of the License at
//
// 	http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
// WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
// License for the specific language governing permissions and limitations
// under the License.

// Package nelmead is s a small implementation of the Nelder-Mead optimization algorithm
// with minimal dependencies. It strictly follows the algorithm as described in:
//
// https://en.wikipedia.org/wiki/Nelder-Mead_method
//
package nelmead

import (
	"sort"
)

// NelderMeadOptimizer is a (tiny) Go implementation of the
// This object contains the basic parameters of the optimizer.
// It is obtained buy using MakeNelderOptimizer.
type NelderMeadOptimizer struct {
	f func([]float64) float64
	// ConvTreshold is the minimimal difference between two successive trials
	// for an optimization to be classified as an improvement.
	ConvTreshold float64
	convCount    int
	// ConvLimit is the number of iterations with no significant improvement in
	// the optimizing function for the solution to be classified as converged.
	ConvLimit int
	Alpha     float64
	Gamma     float64
	Rho       float64
	Sigma     float64
	dim       int
	best      float64
	bestPrev  float64
}

func MakeNelderMeadOptimizer(f func([]float64) float64, dim int) *NelderMeadOptimizer {
	return &NelderMeadOptimizer{
		f:            f,
		ConvTreshold: 10E-15,
		ConvLimit:    12,
		Alpha:        1.0,
		Gamma:        2.0,
		Rho:          0.5,
		Sigma:        0.5,
		dim:          dim,
	}

}

type VertexScore struct {
	X     []float64
	Score float64
}

type SimplexValues []VertexScore

func MakeSimplexValues(start []float64, step float64) SimplexValues {
	dim := len(start)
	spx := make([]VertexScore, dim+1)
	spx[0].X = start

	// Initial values for other vertices, by stepping along each coordinate direction
	for i := 1; i <= dim; i++ {
		v := make([]float64, dim)
		copy(v, start)
		v[i-1] += step
		spx[i].X = v
	}
	return spx
}

// SimplexValues.Centroid calculates the centroid of all vertices of the simplex,
// except the last vertex given by spx[dim]
func (spx SimplexValues) Centroid() []float64 {
	xctr := make([]float64, len(spx)-1)
	for j := 0; j < len(spx)-1; j++ { // iter dim
		for i := 0; i < len(spx)-1; i++ { // iter points
			xctr[j] += spx[i].X[j]
		}
		xctr[j] /= float64(len(spx) - 1)
	}
	return xctr
}

// SimplexValues.Sort sorts the vertices according to the yield function value, from
// lowest value (best score) to highest (worst score)
func (spx SimplexValues) Sort() {
	sort.Slice(spx, func(i, j int) bool { return spx[i].Score < spx[j].Score })
}

// SimplexValues.Set is a helper method to set the vertex values and the scores
// in the Simplex array.
func (spx SimplexValues) Set(i int, x []float64, s float64) {
	if 0 <= i && i < len(spx) {
		spx[i].X = x
		spx[i].Score = s
	}
}

// CheckTerminate is used to check if the it still makes sense to continue the
// iteration process. If successive values for the yield function are getting
// too close, the iteration process will stop as apparently a (local) minimum
// is found.
func (nmo *NelderMeadOptimizer) CheckTerminate() bool {
	if nmo.best < (nmo.bestPrev - nmo.ConvTreshold) {
		nmo.convCount = 0
		nmo.bestPrev = nmo.best
	} else {
		nmo.convCount += 1
	}

	if nmo.convCount >= nmo.ConvLimit {
		return true
	} else {
		return false
	}

}

// MakePoint creates a new point, using a vector between two points. This is a
// helper method, to simplify the generation of the new points in the different
// stages of the optimization process.
func MakePoint(factor float64, x0, xm []float64) []float64 {
	xr := make([]float64, len(x0))
	for i := 0; i < len(x0); i++ {
		xr[i] = x0[i] + factor*(xm[i]-x0[i])
	}
	return xr
}

// Optimize is the key method of this package, and in this method the actual iteration
// process of the Nelder-Mead algorithm is implemented. For more detailed information see the
// Wikimedia entry at: https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method.
func (nmo *NelderMeadOptimizer) Optimize(max_iter int, start []float64, step float64) ([]float64, float64, bool) {
	// initializations
	spx := MakeSimplexValues(start, step)
	for i := 0; i <= nmo.dim; i++ {
		spx[i].Score = nmo.f(spx[i].X)
	}
	nmo.bestPrev = spx[0].Score

	for iters := 0; iters < max_iter; iters++ {
		//
		// (1) order results, best score at lowest index of SimplexValues
		spx.Sort()
		nmo.best = spx[0].Score

		// check termination conditions
		if nmo.CheckTerminate() {
			return spx[0].X, spx[0].Score, true
		}

		// (2) calculate  centroid for dim points, exclude dim+1
		xctr := spx.Centroid()

		// (3) Reflection
		xr := MakePoint(-nmo.Alpha, xctr, spx[nmo.dim].X) // negative alpha to use general MakePoint function
		rscore := nmo.f(xr)

		switch {
		case spx[0].Score <= rscore && rscore < spx[nmo.dim-1].Score:
			spx.Set(nmo.dim, xr, rscore)
		case rscore < spx[0].Score: // step (4) Expansion
			xe := MakePoint(nmo.Gamma, xctr, xr)
			escore := nmo.f(xe)
			if escore < rscore {
				spx.Set(nmo.dim, xe, escore)
			} else {
				spx.Set(nmo.dim, xr, rscore)
			}
		default: // step (5)
			xc := MakePoint(nmo.Rho, xctr, spx[nmo.dim].X)
			cscore := nmo.f(xc)
			if cscore < spx[nmo.dim].Score {
				spx.Set(nmo.dim, xc, cscore)
			} else { // step (6)
				for i := 1; i <= nmo.dim; i++ {
					x := MakePoint(nmo.Sigma, spx[0].X, spx[i].X)
					spx.Set(i, x, nmo.f(spx[i].X))
				}

			}
		}
	}
	return spx[0].X, spx[0].Score, false
}
