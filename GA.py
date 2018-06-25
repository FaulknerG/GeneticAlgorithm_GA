# -*- coding: utf-8 -*-
#chromosome 染色体
#individuals 个人
#elitist 	精英
#mutation 	变异
#interval 	间歇

import math, random
class Population:
	def __init__ (self, size, chrom_size, cp, mp, gen_max):
		#种群信息
		self.individuals = []			# 个体集合 每个个体为[x,y]
		self.fitness = []				# 个体适应度集合
		self.selector_probability = []	# 个体选择概率集合
		self.new_individuals = []		# 新一代个体集合

		self.elitist = {'chromosome':[0,0], 'fitness':0, 'age':0} # 最佳个体的信息，初始化为0

		self.size = size # 种群所包含的个体数
		self.chromosome_size = chrom_size 	# 个体的染色体长度
		self.crossover_probability = cp 	# 个体之间的交叉概率
		self.mutation_probability = mp 		# 个体之间的变异概率

		self.generation_max = gen_max 	#种群进化的最大世代数
		self.age = 0	#种群当前所处世代

		#随机产生初始个体集，并将新一代个体、适应度、选择概率等集合初始化为0
		v = 2 ** self.chromosome_size - 1
		for i in range(self.size):
			self.individuals.append([random.randint(0, v), random.randint(0, v)]) #每个个体由x,y两个元素组成
			self.new_individuals.append([0,0])
			self.fitness.append(0)
			self.selector_probability.append(0)

	#将一个染色体chromosome映射为区间interval之内的数值。
	def decode (self, interval, chromosome):
		d = interval[1] - interval[0]
		n = float (2 ** self.chromosome_size - 1)
		return (interval[0] + chromosome * d / n)

	#适应度函数，可以根据个体的两个染色体计算出该个体的适应度
	#n: sqrt(sin(x2+y2)) -0.5
    #

	def fitness_func (self, chrom1, chrom2):
		interval = [-1.0, 2.0]
		#把两个染色体映射为-10到10之间的数，代入函数计算函数值（适应度）
		(x,y) = (self.decode (interval, chrom1),
				 self.decode (interval, chrom2))
		#n = lambda x, y: math.sin (math.sqrt(x*x + y*y)) ** 2 - 0.5
		#d = lambda x, y: (1 + 0.001 * (x*x + y*y)) ** 2
		#func = lambda x, y: 0.5 - n (x,y) / d (x,y)
		func = lambda x, y: x * math.sin(10 * math.pi * x) + 2
		return func (x, y)		#返回函数值，即适应度



	#评估种群中的个体集合self.individuals 中各个个体的适应度
	#	即将各个个体的两条染色体代入fitness_func函数中
	#结果保存到self.fitness列表中
	#	将个体适应度除以所有个体适应度之和，即得个体的生存概率，同时可以计算出累加概率
	def evaluate(self):
		sp = self.selector_probability		#个体选择概率集合,即生存概率
		#对每个个体计算适应度
		for i in range (self.size):
			self.fitness[i] = self.fitness_func (self.individuals[i][0],
												 self.individuals[i][1])
		ft_sum = sum (self.fitness)
		for i in range (self.size):
			sp[i] = self.fitness[i] / float (ft_sum)
		for i in range (1, self.size):
			sp[i] = sp[i] + sp[i-1]		#进行累加，得累加的个体选择概率

	#轮盘赌博机机制，产生一个随机数，若其刚好小于第 i 个个体累加概率，则表示随机选中的结果为 i
	def select (self):
		(t, i) = (random.random(), 0)
		for p in self.selector_probability:
			if p > t:
				break
			i = i + 1
		return i 

	#染色体交叉模拟
	#cross函数可以将两个染色体进行交叉配对，从而生成2个新染色体
	def cross (self, chrom1, chrom2):
		p = random.random ()
		n = 2 ** self.chromosome_size - 1 #个体染色体长度，n为组合产生的个数
		if chrom1 != chrom2 and p < self.crossover_probability:	
			#小于交叉概率，即可进行交叉
			t = random.randint (1, self.chromosome_size - 1)	#随机选取一个位置
			mask = n << t 
			(r1, r2) = (chrom1 & mask, chrom2 & mask)
			mask = n >> (self.chromosome_size - t)
			(l1, l2) = (chrom1 & mask, chrom2 & mask)
			(chrom1, chrom2) = (r1 + l2, r2 + l1)
		return (chrom1, chrom2)

	#染色体变异模拟
	#mutate将一个染色体按照变异概率进行单点变异：
	def mutate (self, chrom):
		p = random.random ()
		if p < self.mutation_probability:
			t = random.randint (1, self.chromosome_size)
			mask1 = 1 << (t - 1)
			mask2 = chrom & mask1
			if mask2 > 0:
				chrom = chrom & (~mask2)
			else:
				chrom = chrom ^ mask1
		return chrom


	#进化过程
	#evolve实现种群的一代进化计算，分为三个步骤
	#1.使用evaluate评估当前种群的适应度，并计算各个体的选择概率
	#2.对于数量为self.size的self.individuals集合，循环self.size/2次，每次
	#	从self.individuals中选出两个个体，对其进行交叉和变异操作，结果保存到self.new_individuals中
	#3.用种群进化生成的新个体集合self.new_individuals替换当前个体集合
	def evolve (self):
		indvs = self.individuals
		new_indvs = self.new_individuals 

		#计算适应度以及选择概率
		self.evaluate ()

		#进化操作
		i = 0
		while True:
			#选择两个个体，进行交叉与变异，产生新的种群
			idv1 = self.select ()
			idv2 = self.select ()

			#交叉
			(idv1_x, idv1_y) = (indvs[idv1][0], indvs[idv1][1])
			(idv2_x, idv2_y) = (indvs[idv2][0], indvs[idv2][1])
			(idv1_x, idv2_x) = self.cross (idv1_x, idv2_x)
			(idv1_y, idv2_y) = self.cross (idv1_y, idv2_y)

			#变异
			(idv1_x, idv1_y) = (self.mutate (idv1_x), self.mutate (idv1_y))
			(idv2_x, idv2_y) = (self.mutate (idv2_x), self.mutate (idv2_y))

			new_indvs[i][0], new_indvs[i][1] 	  = (idv1_x, idv1_y)
			new_indvs[i+1][0], new_indvs[i+1][1] = (idv2_x, idv2_y)

			#判断进化过程是否结束
			i = i + 2
			if i >= self.size:
				break

		#更新换代
		for i in range (self.size):
			self.individuals[i][0] = self.new_individuals[i][0]
			self.individuals[i][1] = self.new_individuals[i][1]


	#循环调用evolve函数，就可以产生一个种群进化的过程.
	#run函数根据种群的最大化世代数设定了一个循环。
	def run(self):
		for i in range (self.generation_max):
			self.evolve ()
			print (i, max (self.fitness), sum (self.fitness)/self.size, min (self.fitness))

if __name__ == '__main__':
	#种群个体数50，染色体长度25，交叉概率0.8，变异概率0.1，进化最大世代数150
	pop = Population(50, 24, 0.8, 0.1, 150)
	pop.run ()


		
