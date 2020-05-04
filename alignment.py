def distance(a, b): #Левенштейн
    "Calculates the Levenshtein distance between a and b."
    n, m = len(a), len(b)
    if n > m:
        # Make sure n <= m, to use O(min(n, m)) space
        a, b = b, a
        n, m = m, n

    current_row = range(n + 1)  # Keep current and previous row, not entire matrix
    for i in range(1, m + 1):
        previous_row, current_row = current_row, [i] + [0] * n
        for j in range(1, n + 1):
            add, delete, change = previous_row[j] + 1, current_row[j - 1] + 1, previous_row[j - 1]
            if a[j - 1] != b[i - 1]:
                change += 1
            current_row[j] = min(add, delete, change)

    return current_row[n]

def ALIGN_NW (x, y):
    ALIGNx, ALIGNy = [], []
    N, M = len(x), len(y)
    A, path, W = [0] * (N + 1), [0] * (N), [0] * (N) # Матрица штрафов
    
    for i in range(N + 1):
        A[i] = [0] * (M + 1)
    
    for i in range(N):
        path[i], W[i] = [0] * (M), [0] * (M)
    
    for i in range(N + 1):
        for j in range(M + 1):
            A[0][j], A[i][0] = -2 * (j), -2 * (i)
    
    for i in range(1, N + 1):
        for j in range(1, M + 1):
            d = A[i - 1][j - 1]
            if x[i - 1] == y[j - 1]:
                d += 1
            if x[i - 1] != y[j - 1]:
                d -= 1
            v_put = max(A[i][j - 1] - 2, A[i - 1][j] - 2, d)
            A[i][j] += v_put
    W[0][0], path[0][0] = A[1][1], 'd'
    
    for i in range(N):
        for j in range(M):
            
            if (i == 0 and j != 0):
                W[i][j] = W[i][j - 1] + A[i][j]
                path[i][j] = 'l'
            
            if (j == 0 and i != 0):
                path[i][j] = 'u'
                W[i][j] = W[i - 1][j] + A[i][j]
            
            if (i != 0 and j != 0):
                s = max(W[i][j - 1], W[i - 1][j], W[i - 1][j - 1])
                if (s == W[i][j - 1]):
                    path[i][j] = 'l'
                if (s == W[i - 1][j]):
                    path[i][j] = 'u'
                if (s == W[i - 1][j-1]):
                    path[i][j] = 'd'
                W[i][j] = s + A[i][j]
    
    i, j = 0, 0
    ALIGNx.append(x[i])
    ALIGNy.append(y[j])
    
    while (i < N - 1 and j < M - 1):
        if (i < N - 1 and j < M - 1):
            MAX = max(W[i + 1][j + 1], W[i][j + 1], W[i + 1][j])
            
            if MAX == W[i + 1][j + 1]:
                i += 1
                j += 1
                ALIGNx.append(x[i])
                ALIGNy.append(y[j])
                
            else:
                if MAX == W[i][j + 1]:
                    j += 1
                    ALIGNx.append('-')
                    ALIGNy.append(y[j])
                else:
                    if MAX == W[i + 1][j]:
                        i += 1
                        ALIGNx.append(x[i])
                        ALIGNy.append('-')
                        
        if (i == N - 1 and j != M - 1):
            j += 1
            ALIGNx.append('-')
            ALIGNy.append(y[j])
            
        if (j == M - 1 and i != N - 1):
            i += 1
            ALIGNx.append(x[i])
            ALIGNy.append('-')

    x_end, y_end = '', ''
    for i in range(len(ALIGNx)):
        x_end += ALIGNx[i]
    for j in range(len(ALIGNy)):
        y_end += ALIGNy[j]
            
    return distance(x_end, y_end), x_end, y_end

def FASTA(x, y):
    N, M = len(x), len(y)
    A = ['-'] * N
    
    for i in range(N):
        A[i] = [-1] * M
        
    for i in range(N):
        for j in range(M):
            if x[i] == y[j]:
                A[i][j] = 0
                
    for i in range(N):
        for j in range(M):
            
            if (i != 0 and j != 0):
                if (A[i][j] == 0 and (A[i - 1][j - 1] == 0 or A[i - 1][j - 1] == 1)):
                    A[i][j] = 1
                    
            if (i != N - 1 and j != M - 1):
                if(A[i][j] == 0 and (A[i + 1][j + 1] == 0 or A[i + 1][j + 1] == 1)):
                    A[i][j] = 1
                    
    for i in range(1, N):
        for j in range(1, M):
            if (A[i][j] != 0 and A[i-1][j-1] != -1 and A[i][j] != -1):
                A[i][j] += A[i-1][j-1]
                
    index = []
    k, h, start = 0, 0, 0
    
    for i in range(N):
        if len(index) > 0:
            start = index[h-1][1]
        for j in range(start, M):
            if A[i][j] >= 1:
                if start > 0:
                    if j <= index[h-1][1] or j > index[h-1][1]+1 :
                        k += 1
                index.append([i, j, k])
                h += 1
                break

    ALIGNx, ALIGNy = [], []
    i, j = 0, 0

    for k in range(len(index)):
        
        if (i < index[k][0]):
            while(i < index[k][0]):
                ALIGNx.append(x[i])
                ALIGNy.append('-')
                i += 1
        
        if (j < index[k][1]):
            while(j < index[k][1]):
                ALIGNy.append(y[j])
                ALIGNx.append('-')
                j += 1
        
        if (index[k][0] == index[k-1][0]):
            ALIGNx.append('-')
            ALIGNy.append(y[j])
            i += 1
            j += 1
            
        if (index[k][1] == index[k-1][1]):
            ALIGNy.append('-')
            ALIGNx.append(x[i])
            i += 1
            j += 1
        
        if (i == index[k][0] and j == index[k][1] and index[k][0] != index [k - 1][0] and index[k][1] != index[k - 1][1]):
            ALIGNx.append(x[i])
            ALIGNy.append(y[j])
            i += 1
            j += 1

    x_end, y_end = '', ''
    for i in range(len(ALIGNx)):
        x_end += ALIGNx[i]
    for j in range(len(ALIGNy)):
        y_end += ALIGNy[j]
                   
    return distance(x_end, y_end), x_end, y_end

def BLAST (x, y, w):
    ALIGNx, ALIGNy, strings = '', '', []
    N, M, s = len(x), len(y), len(y) % w
    for i in range(M - w + 1):
        pause = ''
        for j in range(5):
            pause += y[i + j]
        strings.append([pause, i, j + i])
    i, j = 0, 0
    for substr in strings:
        
        left_end, right_end = substr[1], substr[2]
        substring_in_y = substr[0]
        if substring_in_y in x and j < left_end or j == left_end == 0:
            
            init = x.find(substr[0])
            final = init + len(substr[0]) - 1
            substring_in_x = substr[0]
            
            while distance(substring_in_x,substring_in_y) <= 2:
                
                flag_x, flag_y = substring_in_x, substring_in_y
                
                if left_end > 0:
                    left_end -= 1
                    substring_in_y = str(y[left_end]) + str(substring_in_y)
                    if init > 0:
                        init -= 1
                        substring_in_x = str(x[init]) + substring_in_x
                    if init == 0:
                        substring_in_x = '-' + substring_in_x
                
                if right_end < M - 1:
                    right_end += 1
                    substring_in_y = str(substring_in_y) + str(y[right_end])
                    if final < N - 1:
                        final += 1
                        substring_in_x += str(x[final])
                    if final == N - 1:
                        substring_in_x += '-'
                
                if right_end == M - 1 and final < N - 1:
                    substring_in_y += '-'
                    final += 1
                    substring_in_x += str(x[final])
                        
                if left_end == 0 and init > 0:
                    init -= 1
                    substring_in_y = '-' + str(substring_in_y)
                    substring_in_x = str(x[init]) + substring_in_x
                        
                if (flag_x == substring_in_x or flag_y == substring_in_y):
                    break
    
            ALIGNx += substring_in_x            
            ALIGNy += substring_in_y
            j = right_end
    if final < N:
        for i in range(final+1, len(x)):
            ALIGNx += str(x[i])
            ALIGNy += '-'    
    return distance(ALIGNx, ALIGNy), ALIGNx, ALIGNy

