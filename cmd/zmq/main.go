package main

import (
	"fmt"
	"log"
	"math/rand"
	"time"

	"github.com/zeromq/goczmq"
)

const (
	timeout = time.Duration(2500) * time.Millisecond
	retries = 3
	serverEndpoint = "tcp://localhost:5555"
)


// func startClient() {
// 	client, err := goczmq.NewReq(serverEndpoint)
// 	if err != nil {
// 		log.Fatal(err)
// 	}
// 	defer client.Destroy()
//
// 	for triesLeft, try := retries, 0; triesLeft > 0; triesLeft-- {
// 		client.SendMessage([]byte(strconv.Itoa(try));
//
// 		for expectReply := true; expectReply {
// 			poller, err := NewPoller(client)
// 			if err != nil {
// 				panic(err)
// 			}
//
// 			sock := poller.wait(timeout)
// 			if sock == nil {
// 				// No response
// 				if triesLeft {
//
// 				} else {
//
// 				}
// 			}
//
//
// 		}
//
// 	}
//
// }

func workerTask(identity string) {
	dealer, err := goczmq.NewDealer("tcp://127.0.0.1:5555")
	if err != nil {
		log.fatal(err)
	}
	defer dealer.destroy()

	dealer.setidentity(identity)

	total := 0
	for {
		err := dealer.SendMessage([][]byte{ []byte(""), []byte("Ready!") })
		if err != nil {
			print(err)
		}

		// TODO how to handle errors here
		workload, _ := dealer.RecvMessage()
		
		for i, part := range workload {
			fmt.Printf("[Dealer (%s)] Part %d: '%s'\n", identity, i, string(part))
		}
		if string(workload[1]) == "Fired!" {
			id := dealer.Identity()
			fmt.Printf("Completed: %d tasks (%s)\n", total, id)
			break
		}

		total += 1
		msec := rand.Intn(1000)
		time.Sleep(time.Duration(msec) * time.Millisecond)
	}
}

func main() {
	router, err := goczmq.NewRouter("tcp://0.0.0.0:5555")
	if err != nil {
		log.Fatal(err)
	}
	defer router.Destroy()

	rand.Seed(time.Now().Unix())

	for i := 0; i < NBR_WORKERS; i++ {
		go workerTask(fmt.Sprintf("Worker #%d", i))
	}

	end_time := time.Now().Unix() + 5
	workers_fired := 0

	for {
		//  Next message gives us least recently used worker
		parts, err := router.RecvMessage()
		if err != nil {
			print(err)
		}

		for i, part := range parts {
			fmt.Printf("Router: Part %d: '%s'\n", i, string(part))
		}
		identity := parts[0]
		now := time.Now().Unix()
		if now < end_time {
			router.SendMessage([][]byte{identity, []byte(""), []byte("Work harder")})
		} else {
			router.SendMessage([][]byte{identity, []byte(""), []byte("Fired!")})
			workers_fired++
			if workers_fired == NBR_WORKERS {
				break
			}
		}
	}
}
