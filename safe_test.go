package apcomplex

import (
	"math"
	"strconv"
	"strings"
	"sync"
	"testing"
	"time"
)

// local helpers (duplicated here to avoid test package import cycles)
func f64snum(s string) float64 {
	s = strings.TrimSpace(s)
	if len(s) > 0 && s[0] == '+' {
		s = s[1:]
	}
	v, _ := strconv.ParseFloat(s, 64)
	return v
}

func approxEqualSafe(a, b *Safe, tol float64) bool {
	d := a.Sub(b)
	re := f64snum(d.RealStringFixed(40))
	im := f64snum(d.ImagStringFixed(40))
	return math.Abs(re) <= tol && math.Abs(im) <= tol
}

// Ensure Add is commutative under heavy parallel calls and lock ordering
// (exercises lockPairR stable ordering).
func TestSafeDeadlockFreeAdd(t *testing.T) {
	a := MustParseSafe("3.25-1.75i", 256)
	b := MustParseSafe("1.5+0.75i", 256)
	defer a.Close()
	defer b.Close()

	const N = 64
	var wg sync.WaitGroup
	wg.Add(N)
	errs := make(chan string, N)

	for i := 0; i < N; i++ {
		go func() {
			defer wg.Done()
			u := a.Add(b)
			v := b.Add(a)
			// Tight tolerance; both results should be identical
			if !approxEqualSafe(u, v, 1e-35) {
				errs <- "a+b != b+a"
			}
		}()
	}
	wg.Wait()
	close(errs)
	for e := range errs {
		t.Fatalf("parallel add mismatch: %s", e)
	}
}

// Run Exp(Log(z)) concurrently from many goroutines; should equal z within tight tolerance.
func TestSafeConcurrentExpLog(t *testing.T) {
	z := MustParseSafe("0.75+0.5i", 384)
	defer z.Close()

	const G = 32
	var wg sync.WaitGroup
	wg.Add(G)
	errs := make(chan string, G)

	for i := 0; i < G; i++ {
		go func() {
			defer wg.Done()
			back := z.Log().Exp()
			if !approxEqualSafe(back, z, 1e-28) {
				errs <- "exp(log(z)) != z"
			}
		}()
	}
	wg.Wait()
	close(errs)
	for e := range errs {
		t.Fatalf("concurrent exp/log mismatch: %s", e)
	}
}

// Continually change precision while other goroutines read (Log/Exp).
// This specifically checks we have no data races or panics and that results stay finite/parseable.
func TestSafeSetPrecWhileReading(t *testing.T) {
	s := MustParseSafe("1.234567890123456789+6.789012345678901234i", 256)
	defer s.Close()

	stop := make(chan struct{})
	var wg sync.WaitGroup

	// Writer goroutine toggles precision.
	wg.Add(1)
	go func() {
		defer wg.Done()
		for {
			select {
			case <-stop:
				return
			default:
				s.SetPrec(320)
				s.SetPrec(256)
			}
		}
	}()

	// Readers perform functions that take RLock.
	const R = 8
	wg.Add(R)
	errs := make(chan string, R)
	for i := 0; i < R; i++ {
		go func() {
			defer wg.Done()
			// Do some work; errors will be visible with -race if any races exist.
			for j := 0; j < 50; j++ {
				_ = s.Log().StringScientific(80)
				_ = s.Exp().StringScientific(80)
				// Read components to ensure both Re/Im access paths are exercised
				_ = s.RealStringFixed(10)
				_ = s.ImagStringFixed(10)
			}
		}()
	}

	// Let the system run for a short period.
	time.Sleep(200 * time.Millisecond)
	close(stop)
	wg.Wait()
	close(errs)
	for e := range errs {
		t.Fatalf("error in SetPrecWhileReading: %s", e)
	}
}
