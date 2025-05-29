import io

data = [[str(k) for k in range(10)] for i in range(2048)]

with io.open('test.txt', 'w', buffering=524_288) as f:
    for line in data:
        f.write(f"{'\t'.join(line)}\n")