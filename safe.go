package apcomplex

import (
	"sync"
	"unsafe"
)

// Safe wraps a *Complex with a mutex so multiple goroutines can operate on it safely.
// All operations return NEW Safe results; the wrapped value is never mutated externally.
type Safe struct {
	mu sync.RWMutex
	c  *Complex
}

// NewSafe allocates a new Safe complex with the given precision in bits.
func NewSafe(bits uint) *Safe { return &Safe{c: New(bits)} }

// WrapSafe wraps an existing *Complex. After wrapping, do NOT use the raw *Complex concurrently.
func WrapSafe(c *Complex) *Safe { return &Safe{c: c} }

// Close releases resources of the underlying Complex.
func (s *Safe) Close() { s.mu.Lock(); s.c.Close(); s.mu.Unlock() }

// Prec reads the precision (bits).
func (s *Safe) Prec() uint { s.mu.RLock(); p := s.c.prec; s.mu.RUnlock(); return p }

// SetPrec updates precision (rounding value).
func (s *Safe) SetPrec(bits uint) { s.mu.Lock(); s.c.SetPrec(bits); s.mu.Unlock() }

// String/format helpers (read-only)
func (s *Safe) StringFixed(d int) string {
	s.mu.RLock()
	out := s.c.StringFixed(d)
	s.mu.RUnlock()
	return out
}
func (s *Safe) StringScientific(d int) string {
	s.mu.RLock()
	out := s.c.StringScientific(d)
	s.mu.RUnlock()
	return out
}
func (s *Safe) RealStringFixed(d int) string {
	s.mu.RLock()
	out := s.c.RealStringFixed(d)
	s.mu.RUnlock()
	return out
}
func (s *Safe) ImagStringFixed(d int) string {
	s.mu.RLock()
	out := s.c.ImagStringFixed(d)
	s.mu.RUnlock()
	return out
}

// Unsafe returns the underlying *Complex. Use with care (no internal locking).
func (s *Safe) Unsafe() *Complex { return s.c }

// lockPairR acquires read locks on a and b in a stable address order to avoid deadlocks.
func lockPairR(a, b *Safe) (unlock func()) {
	if a == b {
		a.mu.RLock()
		return func() { a.mu.RUnlock() }
	}
	ap := uintptr(unsafe.Pointer(a))
	bp := uintptr(unsafe.Pointer(b))
	if ap < bp {
		a.mu.RLock()
		b.mu.RLock()
		return func() { b.mu.RUnlock(); a.mu.RUnlock() }
	}
	b.mu.RLock()
	a.mu.RLock()
	return func() { a.mu.RUnlock(); b.mu.RUnlock() }
}

// --- Non-mutating arithmetic: each returns a NEW Safe result ---

func (a *Safe) Neg() *Safe {
	a.mu.RLock()
	res := NewSafe(a.c.prec)
	res.c.Neg(a.c)
	a.mu.RUnlock()
	return res
}

func (a *Safe) Conj() *Safe {
	a.mu.RLock()
	res := NewSafe(a.c.prec)
	res.c.Conj(a.c)
	a.mu.RUnlock()
	return res
}

func (a *Safe) Inv() *Safe {
	a.mu.RLock()
	res := NewSafe(a.c.prec)
	res.c.Inv(a.c)
	a.mu.RUnlock()
	return res
}

func (a *Safe) Add(b *Safe) *Safe {
	unlock := lockPairR(a, b)
	defer unlock()
	p := a.c.prec
	if b.c.prec > p {
		p = b.c.prec
	}
	res := NewSafe(p)
	res.c.Add(a.c, b.c)
	return res
}

func (a *Safe) Sub(b *Safe) *Safe {
	unlock := lockPairR(a, b)
	defer unlock()
	p := a.c.prec
	if b.c.prec > p {
		p = b.c.prec
	}
	res := NewSafe(p)
	res.c.Sub(a.c, b.c)
	return res
}

func (a *Safe) Mul(b *Safe) *Safe {
	unlock := lockPairR(a, b)
	defer unlock()
	p := a.c.prec
	if b.c.prec > p {
		p = b.c.prec
	}
	res := NewSafe(p)
	res.c.Mul(a.c, b.c)
	return res
}

func (a *Safe) Div(b *Safe) *Safe {
	unlock := lockPairR(a, b)
	defer unlock()
	p := a.c.prec
	if b.c.prec > p {
		p = b.c.prec
	}
	res := NewSafe(p)
	res.c.Div(a.c, b.c)
	return res
}

func (a *Safe) Pow(b *Safe) *Safe {
	unlock := lockPairR(a, b)
	defer unlock()
	p := a.c.prec
	if b.c.prec > p {
		p = b.c.prec
	}
	res := NewSafe(p)
	res.c.Pow(a.c, b.c)
	return res
}

// Elementary / transcendental (read one, produce new)
func (a *Safe) Sqrt() *Safe {
	a.mu.RLock()
	res := NewSafe(a.c.prec)
	res.c.Sqrt(a.c)
	a.mu.RUnlock()
	return res
}
func (a *Safe) Exp() *Safe {
	a.mu.RLock()
	res := NewSafe(a.c.prec)
	res.c.Exp(a.c)
	a.mu.RUnlock()
	return res
}
func (a *Safe) Log() *Safe {
	a.mu.RLock()
	res := NewSafe(a.c.prec)
	res.c.Log(a.c)
	a.mu.RUnlock()
	return res
}
func (a *Safe) Sin() *Safe {
	a.mu.RLock()
	res := NewSafe(a.c.prec)
	res.c.Sin(a.c)
	a.mu.RUnlock()
	return res
}
func (a *Safe) Cos() *Safe {
	a.mu.RLock()
	res := NewSafe(a.c.prec)
	res.c.Cos(a.c)
	a.mu.RUnlock()
	return res
}
func (a *Safe) Tan() *Safe {
	a.mu.RLock()
	res := NewSafe(a.c.prec)
	res.c.Tan(a.c)
	a.mu.RUnlock()
	return res
}
func (a *Safe) Asin() *Safe {
	a.mu.RLock()
	res := NewSafe(a.c.prec)
	res.c.Asin(a.c)
	a.mu.RUnlock()
	return res
}
func (a *Safe) Acos() *Safe {
	a.mu.RLock()
	res := NewSafe(a.c.prec)
	res.c.Acos(a.c)
	a.mu.RUnlock()
	return res
}
func (a *Safe) Atan() *Safe {
	a.mu.RLock()
	res := NewSafe(a.c.prec)
	res.c.Atan(a.c)
	a.mu.RUnlock()
	return res
}
func (a *Safe) Sinh() *Safe {
	a.mu.RLock()
	res := NewSafe(a.c.prec)
	res.c.Sinh(a.c)
	a.mu.RUnlock()
	return res
}
func (a *Safe) Cosh() *Safe {
	a.mu.RLock()
	res := NewSafe(a.c.prec)
	res.c.Cosh(a.c)
	a.mu.RUnlock()
	return res
}
func (a *Safe) Tanh() *Safe {
	a.mu.RLock()
	res := NewSafe(a.c.prec)
	res.c.Tanh(a.c)
	a.mu.RUnlock()
	return res
}
func (a *Safe) Asinh() *Safe {
	a.mu.RLock()
	res := NewSafe(a.c.prec)
	res.c.Asinh(a.c)
	a.mu.RUnlock()
	return res
}
func (a *Safe) Acosh() *Safe {
	a.mu.RLock()
	res := NewSafe(a.c.prec)
	res.c.Acosh(a.c)
	a.mu.RUnlock()
	return res
}
func (a *Safe) Atanh() *Safe {
	a.mu.RLock()
	res := NewSafe(a.c.prec)
	res.c.Atanh(a.c)
	a.mu.RUnlock()
	return res
}

// Constructors from strings
func ParseSafe(s string, prec uint) (*Safe, error) {
	z, err := Parse(s, prec)
	if err != nil {
		return nil, err
	}
	return WrapSafe(z), nil
}

func MustParseSafe(s string, prec uint) *Safe {
	return WrapSafe(MustParse(s, prec))
}
