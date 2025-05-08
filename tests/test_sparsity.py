from marbel.data_generations import add_extra_sparsity, calc_zero_ratio
import pytest
import numpy as np
import pandas as pd


def generate_sparse_df(rows, cols, init_sparsity, seed):
    rng = np.random.default_rng(seed)
    total = rows * cols
    num_zeros = int(total * init_sparsity)
    num_nonzeros = total - num_zeros

    values = np.concatenate([rng.integers(1, 100, num_nonzeros), np.zeros(num_zeros)])
    rng.shuffle(values)
    values = values.reshape((rows, cols))
    col_names = [f"sample_{i}" for i in range(cols)]
    return pd.DataFrame(values, columns=col_names)

@pytest.mark.parametrize("rows", [5, 10])
@pytest.mark.parametrize("cols", [3, 6])
@pytest.mark.parametrize("init_sparsity", [0.0, 0.1, 0.2, 0.3, 0.4])
@pytest.mark.parametrize("target_sparsity", [0.2, 0.3, 0.5, 0.7])
def test_add_extra_sparsity_random_seed(rows, cols, init_sparsity, target_sparsity):
    seed = np.random.randint(0, 1_000_000)
    df = generate_sparse_df(rows, cols, init_sparsity, seed)
    df = df[(df != 0).any(axis=1)]
    init_sparsity_val = calc_zero_ratio(df)
    max_possible_sparsity = (df.shape[0] * df.shape[1] - df.shape[0]) / df.size

    try:
        result = add_extra_sparsity(df.copy(), target_sparsity, seed)
        new_sparsity = calc_zero_ratio(result)

        assert result.shape == df.shape
        assert (
            new_sparsity >= target_sparsity
            or np.isclose(new_sparsity, target_sparsity)
            or new_sparsity > init_sparsity_val
            or np.isclose(new_sparsity, max_possible_sparsity)
        ), (
            f"Seed={seed}, rows={rows}, cols={cols}, "
            f"init={init_sparsity_val:.3f}, target={target_sparsity:.3f}, "
            f"new={new_sparsity:.3f}, max_possible={max_possible_sparsity:.3f}, "
            f"delta={target_sparsity - init_sparsity_val:.3f}"
        )

        assert (result.values >= 0).all()

    except AssertionError:
        print(f"\nTest failed with seed={seed}, init={init_sparsity_val:.3f}, target={target_sparsity:.3f}")
        df.to_csv(f"failed_case_seed_{seed}_rows_{rows}_cols_{cols}.csv", index=False)
        raise
