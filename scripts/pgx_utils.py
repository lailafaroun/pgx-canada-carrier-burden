def weighted_average(values, weights):
    """Calculate a weighted average after validating inputs."""
    if len(values) != len(weights):
        raise ValueError("Values and weights must have the same length.")
    if not values:
        raise ValueError("At least one value is required.")
    if abs(sum(weights) - 1.0) > 1e-6:
        raise ValueError("Weights must sum to 1.")
    if any(not 0 <= value <= 1 for value in values):
        raise ValueError("Values must be between 0 and 1.")
    if any(weight < 0 for weight in weights):
        raise ValueError("Weights cannot be negative.")

    return sum(value * weight for value, weight in zip(values, weights))


def run_tests():
    assert abs(weighted_average([0.1, 0.3], [0.5, 0.5]) - 0.2) < 1e-12
    assert abs(weighted_average([0.2, 0.4], [0.25, 0.75]) - 0.35) < 1e-12
    print("Shared weighted-average tests passed.")


if __name__ == "__main__":
    run_tests()
