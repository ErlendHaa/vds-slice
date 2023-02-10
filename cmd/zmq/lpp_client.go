package main

import (
	"fmt"
	"log"
	"strconv"
	"time"

	"github.com/zeromq/goczmq"
)

const (
	REQUEST_RETRIES = 3
	REQUEST_TIMEOUT = time.Duration(2500) * time.Millisecond
)

func client() {
	socket, err := goczmq.NewReq("tcp://127.0.0.1:5559")
	if err != nil {
		log.Fatal(err)
	}
	defer socket.Destroy()

	log.Println("I: Connecting to server...")
	
	for sequence, retriesLeft := 1, REQUEST_RETRIES; retriesLeft > 0; sequence++ {
		fmt.Printf("I: Sending (%d)\n", sequence)
		socket.SendFrame([]byte(strconv.Itoa(sequence)), 0)

		for expectReply := true; expectReply; {
			poller, err := goczmq.NewPoller(socket)
			if err != nil {
				panic(err) //  Interrupted
			}

			sock := poller.Wait(int(REQUEST_TIMEOUT.Milliseconds()))
			if sock != nil {
				reply, err := sock.RecvMessage();

				if err != nil {
					panic(err) //  Interrupted
				}

				if replyInt, err := strconv.Atoi(string(reply[0])); replyInt == sequence && err == nil {
					fmt.Printf("I: Server replied OK (%d)\n", replyInt)
					retriesLeft = REQUEST_RETRIES
					expectReply = false
				} else {
					fmt.Printf("E: Malformed reply from server: %s", reply)
				}

			} else if retriesLeft--; retriesLeft == 0 {
				fmt.Println("E: Server seems to be offline, abandoning")
				socket.Destroy()
				break
			} else {
				// No response. Close and reopen a new socket and resend message
				socket.Destroy()
				socket, err := goczmq.NewReq("tcp://127.0.0.1:5559")
				if err != nil {
					log.Fatal(err)
				}
				defer socket.Destroy()
				
				fmt.Printf("W: No response from serverm, retrying... (%d)\n", sequence)
				//  Send request again, on new socket
				socket.SendFrame([]byte(strconv.Itoa(sequence)), 0)

			}
		}
	}
}

func main() {
	client()
}
