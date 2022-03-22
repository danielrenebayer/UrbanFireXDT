/*
 * output.h
 *
 * This contains all functions and classes for buffering outputs
 * and writing them to the disk finally.
 *
 */


namespace output {
    void initializeSubstationOutput(int scenario_id);
}


class OutputQueue {
    /*
     * This class represents a queue of data, that should
     * be written to a file.
     */
    public:
        OutputQueue(const char* filepath, unsigned int max_size); // constructs the object and opens the file
        ~OutputQueue(); // destructs the element and closes the file
        void add_to_queue(); // add element for outputting - blocks, if queue is full
        // run thread
            // checks, if elements are in the queue
            // if not: sleep for some time
            // if yes: empty queue and write data to file
    private:
        const unsigned int max_size;
        // FILE file to write
        // QUEUE ...
};
