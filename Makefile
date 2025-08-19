all: build

build: apcomplex.go cmd/main/main.go
	go build ./...

run: apcomplex.go cmd/main/main.go
	go run ./cmd/main

run2: apcomplex.go cmd/main/main.go
	go run ./cmd/main -base '2+2i' -exp '50-50i' -digits -1 -out fixed -prec 1024

test: apcomplex.go cmd/main/main.go apcomplex_test.go
	go test -v

clean:
	go clean -cache -v

fmt:
	go fmt *.go
