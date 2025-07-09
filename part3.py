

def knapsack_0_1(values, weights, capacity):
    n = len(values)
    
    dp = [[0 for _ in range(capacity + 1)] for _ in range(n + 1)]
    
    for i in range(1, n + 1):
        for w in range(1, capacity + 1):
            if weights[i - 1] <= w:
                dp[i][w] = max(dp[i - 1][w], dp[i - 1][w - weights[i - 1]] + values[i - 1])
            else:
                dp[i][w] = dp[i - 1][w]
    
    return dp

values = [5, 4, 7, 7]  
weights = [5, 6, 8, 4]  
capacity = 13           

dp_matrix = knapsack_0_1(values, weights, capacity)

for row in dp_matrix:
    print(row)
