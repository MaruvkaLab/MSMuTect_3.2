import unittest
from src.FullAnalyses.SingleFileBatches import *


class TestSingleFileBatches(unittest.TestCase):
    @staticmethod
    def validate_chunks(chunks: List[Chunk]):
        if len(chunks) == 1:
            return True
        prev_chunk = chunks[0]
        for i in range(1, len(chunks)):
            current_chunk = chunks[i]
            if prev_chunk.end != current_chunk.start-1 or abs(prev_chunk.end + current_chunk.start - prev_chunk.start - current_chunk.end) > 1:
                return False
            prev_chunk = current_chunk
        return True

    def test_chunking(self):
        self.assertTrue(self.validate_chunks(get_chunks(1, 2, 10)))
        self.assertTrue(self.validate_chunks(get_chunks(3, 10, 100)))
        self.assertTrue(self.validate_chunks(get_chunks(2, 0, 1000)))


if __name__ == '__main__':
    unittest.main()
