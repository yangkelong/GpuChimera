import pickle


def saveData(data, file="./data/arr.data"):
    with open(file, "wb") as f:
        pickle.dump(data, f)


def loadData(file="./data/arr.data"):
    with open(file, "rb") as f:
        data = pickle.load(f)
    return data


if __name__ == "__main__":
    print(loadData(r'F:\pythonProject\Output\f6_medium\data18000.fluid'))

