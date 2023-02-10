package main

import (
	"log"

	"github.com/zeromq/goczmq"
)

func queue() {
	proxy := goczmq.NewProxy()
	defer proxy.Destroy()

	if err := proxy.SetBackend(goczmq.Router, "tcp://127.0.0.1:5557"); err != nil {
		panic(err)
	}
	
	if err := proxy.SetFrontend(goczmq.Router, "tcp://127.0.0.1:5559"); err != nil {
		panic(err)
	}

	if err := proxy.Verbose(); err != nil {
		panic(err)
	}
	
	if err := proxy.Resume(); err != nil {
		panic(err)
	}
	
	run := true
	for run { }
	log.Println("exiting...")
}

func main() {
	queue()
}
