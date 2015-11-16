/*

ProducerConsumerQueue.tpp - This file is part of ProducerConsumerQueue (v0.1)

The MIT License (MIT)

Copyright (c) 2014 Lasse Maretty and Jonas Andreas Sibbesen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*/

template <typename Data>
ProducerConsumerQueue<Data>::ProducerConsumerQueue(uint max_buffer_size_in) {

	max_buffer_size = max_buffer_size_in;
	pushed_last = false;
}


template<typename Data> 
void ProducerConsumerQueue<Data>::push(Data data) {

    std::unique_lock<std::mutex> queue_lock(queue_mutex);                              
                
	while (queue.size() == max_buffer_size) {

        producer_cv.wait(queue_lock);
   	}

    queue.push_back(data);
    queue_lock.unlock();

    consumer_cv.notify_one();
}

template<typename Data>
void ProducerConsumerQueue<Data>::pushedLast() {

	std::unique_lock<std::mutex> queue_lock(queue_mutex);
    pushed_last = true;
    queue_lock.unlock();

    consumer_cv.notify_all();
}

template<typename Data>
bool ProducerConsumerQueue<Data>::pop(Data * data) {

    std::unique_lock<std::mutex> queue_lock(queue_mutex);

    while (queue.size() == 0) {

        if (pushed_last) {

            return false;
        }

        consumer_cv.wait(queue_lock);
    }
    
    (*data) = queue.front();
    queue.pop_front();
	queue_lock.unlock();

    producer_cv.notify_one();

    return true;
}