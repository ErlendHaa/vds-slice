package main

import (
	"fmt"
	"log"
	"math/rand"

	"github.com/zeromq/goczmq"
)

const (
	HEARTBEAT_LIVENESS = 3
	HEARTBEAT_INTERVAL = time.Second
)

type Worker struct {
	address []byte
	expiry  time.Time
}

func NewWorker(address string) Worker {
	return Worker {
		address: address,
		expiry:  time.Now().Add(HEARTBEAT_LIVENESS * HEARTBEAT_INTERVAL),
	}
}

type WorkerQueue struct {
	queue *list.List
}


func NewWorkerQueue() WorkerQueue {
	return WorkerQueue {
		queue: list.New()
	}
}

func (w WorkerQueue) Len() int {
	return w.queue.Len()
}

func (w WorkerQueue) Next() Worker {
	element := w.queue.Back()
	worker, _ := element.Value.(Worker)
	w.queue.Remove(element)
	return worker
}

func (w WorkerQueue) Add(worker Worker) {
	for x := w.queue.Front; x != nil; x = x.Next() {
		w := x.Value.(Worker)
		if w.address == worker.address {
			w.queue.Remove(w)
			break
		}
	}

	w.queue.PushBack(w)
}

func (w *WorkerQueue) Purge() {
	now := time.Now()
	for elem := w.queue.Front(); elem != nil; elem = w.queue.Front() {
		if w, _ := elem.Value.(*PPWorker); w.expiry.After(now) {
			break
		}
		w.queue.Remove(elem)
	}
}

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
