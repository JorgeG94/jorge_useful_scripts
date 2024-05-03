from fractions import Fraction

class Resource_Request:
    def __init__(self, queue, ncpus, ngpus, mem, jobfs, minutes = 120):
        self.queue = queue
        self.ncpus = ncpus
        self.ngpus = ngpus 
        self.mem = mem 
        self.jobfs = jobfs
        self.minutes = minutes

class Queue:
    def __init__(self, name, cpus_per_node, gpus_per_node, mem_per_node, jobfs_per_node):
        self.name = name
        self.cpus_per_node = cpus_per_node
        self.gpus_per_node = gpus_per_node
        self.mem_per_node = mem_per_node
        self.jobfs_per_node = jobfs_per_node

    def resource_request(self, ncpus = None, ngpus = None, minutes=120):
        if ncpus == None and ngpus == None:
            return None

        if ncpus and ngpus:
            raise Exception("Requested both cpus and gpus")

        if ncpus and ncpus > self.cpus_per_node:
            raise Exception(f"Requested too many CPUs: queue {self.name}, {ncpus} out of {self.cpus_per_node}")

        if ngpus and ngpus > self.gpus_per_node:
            raise Exception(f"Requested too many CPUs: queue {self.name}, {ngpus} out of {self.gpus_per_node}")

        if ngpus:
            fraction = Fraction(ngpus,self.gpus_per_node)
        elif ncpus:
            fraction = Fraction(ncpus,self.cpus_per_node)

            if self.gpus_per_node > 0 and (self.gpus_per_node*fraction).denominator != 1:
                raise Exception(f"For {self.name} queue, ncpus ({ncpus}) must be a multiple of {self.cpus_per_node//self.gpus_per_node}")


        return Resource_Request(self.name,
                                round(self.cpus_per_node*fraction), 
                                round(self.gpus_per_node*fraction), 
                                round(self.mem_per_node*fraction), 
                                round(self.jobfs_per_node*fraction),
                                minutes)

queues = {}
queues["normal"] = Queue("normal", 48, 0, 190, 400)
queues["normalsr"] = Queue("normalsr", 104, 0, 500, 400)
queues["gpuvolta"] = Queue("gpuvolta", 48, 4, 382, 400)
queues["dgxa100"] = Queue("dgxa100", 128, 8, 2000, 28000)

queues["dgxa100"].resource_request(ncpus=16)
