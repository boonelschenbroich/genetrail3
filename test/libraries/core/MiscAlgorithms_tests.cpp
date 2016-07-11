#include <genetrail2/core/misc_algorithms.h>

#include <gtest/gtest.h>

#include <random>

using namespace GeneTrail;

class MiscAlgorithmsTest : public ::testing::Test
{
};

TEST_F(MiscAlgorithmsTest, testPartialShuffle)
{
	std::mt19937_64 rng(std::random_device{}());

	std::vector<int> data{1, 2, 3, 4, 5, 6};

	// Stress test the algorithm to make sure, that no entries get lost.
	for(int i = 0; i < 1000; ++i) {
		// Determine a random split point
		auto until = rng() % data.size();

		// Shuffle the range between begin and until
		partial_shuffle(data.begin(), data.begin() + until, data.end(), rng);

		// Make sure nothing was lost
		std::sort(data.begin(), data.end());
		ASSERT_EQ(1, data[0]);
		ASSERT_EQ(2, data[1]);
		ASSERT_EQ(3, data[2]);
		ASSERT_EQ(4, data[3]);
		ASSERT_EQ(5, data[4]);
		ASSERT_EQ(6, data[5]);
	}
}

TEST_F(MiscAlgorithmsTest, testSortPermutation)
{
	std::vector<unsigned int> data { 3, 2, 5, 1, 8, 7};

	auto permutation = sort_permutation(data.begin(), data.end(), std::less<size_t>());

	EXPECT_EQ(3u, permutation[0]);
	EXPECT_EQ(1u, permutation[1]);
	EXPECT_EQ(0u, permutation[2]);
	EXPECT_EQ(2u, permutation[3]);
	EXPECT_EQ(5u, permutation[4]);
	EXPECT_EQ(4u, permutation[5]);
}

TEST_F(MiscAlgorithmsTest, testInvertPermutation)
{
	std::vector<size_t> perm { 0, 3, 4, 6, 7, 1, 8, 5, 9, 2 };

	auto inv_perm = invert_permutation(perm);

	EXPECT_EQ(0u, inv_perm[0]);
	EXPECT_EQ(5u, inv_perm[1]);
	EXPECT_EQ(9u, inv_perm[2]);
	EXPECT_EQ(1u, inv_perm[3]);
	EXPECT_EQ(2u, inv_perm[4]);
	EXPECT_EQ(7u, inv_perm[5]);
	EXPECT_EQ(3u, inv_perm[6]);
	EXPECT_EQ(4u, inv_perm[7]);
	EXPECT_EQ(6u, inv_perm[8]);
	EXPECT_EQ(8u, inv_perm[9]);
}

TEST_F(MiscAlgorithmsTest, stressTestInvertPermutation)
{
	std::mt19937_64 rng(std::random_device{}());

	for(unsigned int i = 0; i < 1000; ++i) {
		std::vector<uint64_t> scores(rng() % 1000);

		std::iota(scores.begin(), scores.end(), static_cast<uint64_t>(0));
		std::shuffle(scores.begin(), scores.end(), rng);

		auto perm = sort_permutation(scores.begin(), scores.end(), std::less<uint64_t>());
		auto inv_perm = invert_permutation(perm);

		for(size_t i = 0; i < perm.size(); ++i) {
			EXPECT_EQ(i, perm[inv_perm[i]]);
			EXPECT_EQ(i, inv_perm[perm[i]]);
		}
	}
}
