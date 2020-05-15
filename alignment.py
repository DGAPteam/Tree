def readseq(Sequences):
#_input = open('example/input.fas', 'r')
#Sequences = _input.read()
#_input.close()

    Seq, names_seq = [], []
    Seq_inner = []
    Seq = Sequences.split('>')
    for i in Seq:
        Seq_inner.append(i.split('\n'))
    Seq_inner.remove(Seq_inner[0])
    sequences = []
    for i in Seq_inner:
        s = ''
        names_seq.append(i[0])
        for j in range(1,len(i)):
            s += i[j]
        sequences.append(s)
    if len(sequences[0]) > len(sequences[1]):
        name_seq = names_seq[1]
        name_seqi = names_seq[0]
        seq = sequences[1]
        seqi = sequences[0]
    else:
        name_seq = names_seq[0]
        name_seqi = names_seq[1]
        seq = sequences[0]
        seqi = sequences[1]
    return(seq, seqi)

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

def hamming_distance_nw(s1, s2):
    return sum(ch1 != ch2 and ch1 != '-' and ch2 != '-' for ch1,ch2 in zip(s1,s2))

def hamming_distance_fasta(s1, s2):
    return sum(ch1 != ch2 for ch1,ch2 in zip(s1,s2))

def hamming_distance_blast(s1, s2):
    return sum(ch1 != ch2 for ch1,ch2 in zip(s1,s2))

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
    
    while (i < N - 1 or j < M - 1):
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

    length = max(len(ALIGNx),len(ALIGNy))
    for i in range(length):
        x_end += ALIGNx[i]
    for j in range(length):
        y_end += ALIGNy[j]
    power = max(N, M)

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
            if (A[i][j] == 1 and A[i - 1][j - 1] != -1):
                A[i][j] += A[i - 1][j - 1]
                
    for i in range(N-2, -1, -1):
        for j in range(M-2, -1, -1):
            if (A[i][j] >= 1 and A[i + 1][j + 1] >= 1):
                A[i][j] = A[i + 1][j + 1]
                
    index = []
    h, start_x, start_y = 0, 0, 0
    
    for i in range(start_x, N):
        for j in range(start_y, M):
            if A[i][j] >= 5:
                index.append([i, j])
                h += 1
                start_y = index[h - 1][1] + 1
                start_x = index[h - 1][0] + 1
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

    length = max(len(ALIGNx),len(ALIGNy))
    for i in range(length):
        x_end += ALIGNx[i]
    for j in range(length):
        y_end += ALIGNy[j]
    power = max(N, M)
    while(len(x_end) < power):
        x_end += '-'
    while(len(y_end) < power):
        y_end += '-'
    k, l = 0, 0
    for i in x_end:
        if i != '-':
            k += 1
    for j in y_end:
        if j != '-':
            l += 1
    k += 1
    l += 1
    while k < N:
        x_end += str(x[k])
        k += 1
    while l < M:
        y_end += str(y[l])
        l += 1
    if (len(x_end) < len(y_end)):
        while(len(x_end) < len(y_end)):
            x_end += '-'
    if (len(y_end) < len(x_end)):
        while(len(y_end) < len(x_end)):
            y_end += '-'
    return distance(x_end, y_end), x_end, y_end, k, l

def BLAST (x, y, w):
    ALIGNx, ALIGNy, strings = '', '', []
    N, M, s = len(x), len(y), len(y) % w
    for i in range(M - w + 1):
        pause = ''
        for j in range(w):
            pause += y[i + j]
        strings.append([pause, i, j + i])
    i, j = 0, 0 # в строке x, в строке y
    for substr in strings:
        
        left_end, right_end = substr[1], substr[2] # в строке y
        substring_in_y = substr[0]
        if x.find(substr[0]) > -1 and i <= x.find(substr[0]) and j <= left_end:
            
            init = x.find(substr[0])
            final = init + len(substr[0]) - 1 # init и final в x
            substring_in_x = substr[0]
            
            while hamming_distance_blast(substring_in_x, substring_in_y) <= 1:
                
                flag_x, flag_y = substring_in_x, substring_in_y
                
                if left_end > j:
                    left_end -= 1
                    substring_in_y = str(y[left_end]) + str(substring_in_y)
                    if init > i:
                        init -= 1
                        substring_in_x = str(x[init]) + substring_in_x
                    if init == i:
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
                        
                if left_end == j and init > i:
                    init -= 1
                    substring_in_y = '-' + str(substring_in_y)
                    substring_in_x = str(x[init]) + substring_in_x
                        
                if (flag_x == substring_in_x or flag_y == substring_in_y):
                    break

            if (j < left_end):
                while (j < left_end):
                    ALIGNy += y[j]
                    j += 1
            if (i < init):
                while (i < init):
                    ALIGNx += x[i]
                    i += 1
            
            ALIGNx += substring_in_x            
            ALIGNy += substring_in_y
            j = right_end + 1
            i = final + 1

    x_end, y_end = '', ''

    for i in ALIGNx:
        x_end += i
    for j in ALIGNy:
        y_end += j
        
    k, l = 0, 0
    for i in x_end:
        if i != '-':
            k += 1
    for j in y_end:
        if j != '-':
            l += 1
    
    if k < N:
        while k < N:
            x_end += str(x[k])
            k += 1
    if l < M:
        while l < M:
            y_end += str(y[l])
            l += 1
            
    if (len(x_end) < len(y_end)):
        while(len(x_end) < len(y_end)):
            x_end += '-'
    if (len(y_end) < len(x_end)):
        while(len(y_end) < len(x_end)):
            y_end += '-'
    return distance(x_end, y_end), x_end, y_end, k, l

def GenSearch(x):
    start, end = [], []
    i, j = 0, 0
    while 1:
        l = x.find('ATG', i, len(x) - 1)
        if l == -1:
            break
        if l % 3 == 0:
            start.append(l)
        i = l + 1
    while 1:
        l = x.find('TAA', j, len(x) - 1)
        if l != -1 and l % 3 == 0:
            end.append(l)
        k = x.find('TAG', j, len(x) - 1)
        if k != -1 and k % 3 == 0:
            end.append(k)
        p = x.find('TGA', j, len(x) - 1)
        if p != -1 and p % 3 == 0:
            end.append(p)
        j = max(p, l, k) + 1
        if ( k == l == p == -1):
            break
    genes = []
    k = 0
    for i in start:
        f = 0
        if int(i) > k :
            for j in end:
                if int(j) > i :
                    f = 1
                    s = x[int(i): int(j)]
                    k = int(j)
                    break
        if f == 1:
            genes.append(s)
    return genes
            
def NW_genes(seq, seqi, output_name):
    #output_name = 'output.Needleman-Wunch(Genes)'
    Out = open(output_name, 'w')
    genes1, genes2 = GenSearch(seq), GenSearch(seqi)
    Out.write(name_seq + '\n' + name_seqi + '\n')
    for i in genes1:
        for j in genes2:
            #Для примера
            s1, s2 = i, j
            if (len(s1) > len(s2)):
                while(len(s1) > len(s2)):
                    s2 += '-'
            if (len(s2) > len(s1)):
                while(len(s2) > len(s1)):
                    s1 += '-'
            _distance = distance(s1, s2)
            align = ALIGN_NW(j, i) if len(j) < len(i) else ALIGN_NW(i, j)
            Out.write("длины строк " + str(len(i)) + ' ' + str(len(j)) + ' ' + str(len(align[1])) + ' ' + str(len(align[2])))
            Out.write('\n' + str(_distance/len(s1)) + '\n' + str(hamming_distance_nw(s1,s2)/len(s1)) + '\n' + '\n')
            Out.write(str(hamming_distance_nw(align[1],align[2])/len(align[1])) + '\n' + str(align[0]/len(align[1])) + '\n' + str(align[1]) + '\n' + str(align[2]) + '\n')
    Out.close()
