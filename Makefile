all: build

build: fmt apcomplex.go cmd/main/main.go
	go build ./...

run_main: fmt apcomplex.go cmd/main/main.go
	go run ./cmd/main

run_main2: fmt apcomplex.go cmd/main/main.go
	go run ./cmd/main -base '2+2i' -exp '50-50i' -digits -1 -out fixed -prec 1024

tetrate: fmt apcomplex.go cmd/main/main.go
	go run ./cmd/tetrate 2 5 100

test: fmt apcomplex.go apcomplex_test.go
	go test -v -race ./...

clean:
	go clean -cache -v

fmt:
	go fmt *.go
	go fmt cmd/main/main.go
