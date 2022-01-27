#!/usr/bin/env python3

# based on the chromsweep algorithm from https://github.com/arq5x/bedtools/blob/master/src/utils/chromsweep/chromsweep.cpp

class OverlapFinder:
    def __init__(self, left_file_path, right_file_path, output_path="output.csv", threshold=0.5):
        self.threshold = threshold/100
        self.right_file_path = right_file_path
        self.left_file_path = left_file_path
        self.output_path = output_path

        self.right_fp = open(self.right_file_path, "r")
        self.left_fp = open(self.left_file_path, "r")
        self.output_fp = open(output_path, "w")
        self.output = []

    @staticmethod
    def find_next(iter_):
        try:
            line = next(iter_)
            line = line.strip().split("\t", 3)
            return [line[0], int(line[1]), int(line[2])]
        except:
            return None

    def is_overlapping(self, left_, right_):
        left_len = left_[2] - left_[1]
        overlap = min(right_[2], left_[2]) - max(left_[1], right_[1])
        if left_len == 0:
            return overlap >= 0
        if overlap / left_len >= self.threshold:
            return True
        return False

    def find_overlaps(self):
        ans = []
        right_buffer = []

        # find the first elements of left and right
        current_left = self.find_next(self.left_fp)
        current_right = self.find_next(self.right_fp)

        while current_left is not None:
            # check if chromosomes of left and right match. if they dont match iterate the other till they match
            current_left, current_right, right_buffer, ans = self.sync(current_left, current_right, right_buffer, ans)
            # scan the buffer to find possible candidates for intersection

            right_buffer = self.find_candidates(current_left, right_buffer, ans)

            while current_right is not None and current_right[0] == current_left[0] and not current_right[1] >= current_left[2]:
                right_buffer.append(current_right)
                if self.is_overlapping(current_left, current_right):
                    ans.append(current_right)
                current_right = self.find_next(self.right_fp)

            self.write_ans(current_left, ans)
            ans = []
            current_left = self.find_next(self.left_fp)

        self.write_to_file()

    def find_candidates(self, current_left, right_buffer, ans):
        if current_left is None:
            return right_buffer
        buff = []
        for curr_right in right_buffer:
            if curr_right[0] == current_left[0] and not current_left[1] >= curr_right[2]:
                buff.append(curr_right)
                if self.is_overlapping(current_left, curr_right):
                    ans.append(curr_right)
        return buff

    def sync(self, current_left, current_right, right_buffer, ans):
        if current_right is None or current_left[0] == current_right[0]:
            return current_left, current_right, right_buffer, ans
        if current_left[0] > current_right[0]:
            # left is greater than right
            temp_right = current_right
            while temp_right is not None and temp_right[0] < current_left[0]:
                temp_right = self.find_next(self.right_fp)
            return current_left, temp_right, [], ans
        elif current_left[0] < current_right[0]:
            # right is greater than left
            temp_left = current_left
            while temp_left is not None and temp_left[0] == current_left[0]:
                right_buffer = self.find_candidates(temp_left, right_buffer, ans)
                self.write_ans(temp_left, ans)
                temp_left = self.find_next(self.left_fp)
                ans = []
            while temp_left is not None and temp_left[0] < current_right[0]:
                self.write_ans(temp_left, ans)
                temp_left = self.find_next(self.left_fp)
            return temp_left, current_right, [], ans

    def write_ans(self, left_, ans):
        for right_ in ans:
            assert left_[0] == right_[0], f"Chromosome of left and right dont match"
        if args.j == True:
            for right_ in ans:
                self.output.append(f"{left_[0]}\t{left_[1]}\t{left_[2]}\t{right_[0]}\t{right_[1]}\t{right_[2]}\n")
        else:
            for right_ in ans:
                self.output.append(f"{right_[0]}\t{right_[1]}\t{right_[2]}\n")

    def write_to_file(self):
        self.output_fp.writelines(self.output)

        self.output_fp.close()
        self.right_fp.close()
        self.left_fp.close()
if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i1', help="Input first file", required=True, type=str)
    parser.add_argument('-i2', help="Input second file", required=True, type=str)
    parser.add_argument('-m', help="Threshold value", type=int)
    parser.add_argument("-j", help="Join flag", action='store_true')
    parser.add_argument("-o", help="Name of output file", type=str)
    args = parser.parse_args()

    finder = OverlapFinder(args.i1, args.i2, args.o, args.m)
    finder.find_overlaps()
