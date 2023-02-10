package main

import (
	"fmt"
	"log"
	"math/rand"

	"github.com/zeromq/goczmq"
)

func randomIdentity() string {
	return fmt.Sprintf("%04X-%04X", rand.Intn(0x10000), rand.Intn(0x10000))
}

func queue() {
	frontend, err := goczmq.NewRouter("tcp://127.0.0.1:5559")
	if err != nil {
		log.Fatal(err)
	}
	defer frontend.Destroy()

	backend, err := goczmq.NewRouter("tcp://127.0.0.1:5557")
	if err != nil {
		log.Fatal(err)
	}
	defer backend.Destroy()
	
	frontend.SetIdentity(randomIdentity())
	backend.SetIdentity(randomIdentity())
	
	log.Println(backend.Identity())
	log.Println(frontend.Identity())

	workers := make([][]byte, 0, 0)

	log.Println("I: Queue ready...")
	for {

		poller, err := goczmq.NewPoller(backend, frontend)
		if err != nil {
			panic(fmt.Sprintf("E: could not create poller, got: %s", err)) //Interrupted
		}

		socket := poller.Wait(1)
		if socket == nil {
			continue
		}

		if socket.Identity() == backend.Identity() {
			msg, err := backend.RecvMessage(); 
			if err != nil {
				panic(fmt.Sprintf("E: could not receive on backend, got: %s", err)) //Interrupted
			}

			address := msg[0]
			workers = append(workers, address)

			if reply := msg[2:]; string(reply[0]) != "READY" {
				log.Println("I: routing reply")
				frontend.SendMessage(reply)
			}

			log.Printf("I: Number of workers availible (%d)", len(workers))
		}

		if len(workers) > 0 && socket.Identity() == frontend.Identity() {
			msg, err := frontend.RecvMessage(); 
			if err != nil {
				panic(fmt.Sprintf("E: could not receive on frontend, got: %s", err)) //Interrupted
			}
			
			worker := workers[0]
			workers = workers[1:]
			log.Println("I: routing requst")

			// log.Println("I: routing request (%v)", msg)
			request := append([][]byte{worker, nil}, msg...)
			backend.SendMessage(request)
			log.Printf("I: Number of workers availible (%d)", len(workers))
		}
	}
}

func main() {
	queue()
}
