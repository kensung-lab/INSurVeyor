import random

class RandomPositionGenerator:

    NO_SAMPLING_PADDING = 1000

    def __init__(self, reference_fa, seed, sampling_regions_fname = None):
        random.seed(seed)

        sampling_regions = []
        if sampling_regions_fname:
            with open(sampling_regions_fname) as sr_inf:
                for line in sr_inf:
                    sl = line.split()
                    if len(sl) == 1:
                        chr, start, end = sl[0], self.NO_SAMPLING_PADDING, len(reference_fa[sl[0]])-self.NO_SAMPLING_PADDING
                    else:
                        chr, start, end = sl[0], int(sl[1]), int(sl[2])
                    sampling_regions += [(chr, start, end)]
        else:
            sampling_regions = [(k, self.NO_SAMPLING_PADDING, len(reference_fa[k])-self.NO_SAMPLING_PADDING) for k in reference_fa.keys()]

        self.sampling_regions = sampling_regions
        self.reference_len = sum([r[2]-r[1] for r in sampling_regions])

    def next(self):
        # according to documentation this is supposed to generate within [a, b), but sometimes it generates b for me
        r = random.randint(0, self.reference_len-1)
        for region in self.sampling_regions:
            k, start, stop = region
            if r < stop-start:
                chr = k
                break
            r -= stop-start
        return (chr, r+start)
