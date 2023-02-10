package main

import (
	"fmt"
	"log"
	"math/rand"
	"time"

	"github.com/zeromq/goczmq"
)

func worker() {
	// Random string
	identity := fmt.Sprintf("%04X-%04X", rand.Intn(0x10000), rand.Intn(0x10000))

	socket, err := goczmq.NewReq("tcp://127.0.0.1:5557")
	if err != nil {
		log.Fatal(err)
	}
	defer socket.Destroy()

	socket.SetHeartbeatIvl(1)
	socket.SetIdentity(identity)
	
	rand.Seed(time.Now().Unix())

	//  Tell broker we're ready for work
	fmt.Printf("I: (%s) worker ready\n", identity)
	err = socket.SendFrame([]byte("READY"), 0)
	if err != nil {
		print(err)
	}

	for cycles := 1; ; cycles++ {
		msg, err := socket.RecvMessage();
		if err != nil {
			panic(err)
		}
		
		if cycles > 3 {
			switch r := rand.Intn(5); r {
			case 0:
				fmt.Printf("I: (%s) simulating a crash\n", identity)
				return
			case 1:
				fmt.Printf("I: (%s) simulating CPU overload\n", identity)
				time.Sleep(3 * time.Second)
			}
		}

		fmt.Printf("I: (%s) normal reply (%s)\n", identity, string(msg[0]))
		time.Sleep(1 * time.Second) //  Do some heavy work
		socket.SendMessage(msg)
	}
}

func main() {
	worker()
}
