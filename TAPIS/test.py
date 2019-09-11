dic = {'Tim': 'Scruffy', 'Peter': 'Furry', 'Sally': 'Fluffy'}

words = ["happy", "Tim", "Peter"]
for word in words:
	if dic.get(word):
		print dic[word]