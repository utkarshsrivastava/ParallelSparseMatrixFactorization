 SMF 
 
  Sparse Matrix Factorization (SMF) is a key component in many machine learning problems

and there exist a verity a applications in real-world problems such as recommendation systems,

estimating missing values, gene expression modeling, intelligent tutoring systems (ITSs), etc.

There are different approaches to tackle with SMF rooted in linear algebra and probability

theory. In this project, given an incomplete binary matrix of students’ performances over a

set of questions, estimating the probability of success or fail over unanswered questions is of

interest. This problem is formulated using Maximum Likelihood Estimation (MLE) which leads

to a biconvex optimization problem (this formulation is based on SPARFA [4]). The resulting

optimization problem is a hard problem to deal with due to the existence of many local minima.

On the other hand, when the size of the matrix of students’ performances increase, the existing

algorithms are not successful; therefore, an efficient algorithm is required to solve this problem

for large matrices. In this project, a parallel algorithm (i.e., a parallel version of SPARFA) is

developed to solve the biconvex optimization problem and tested via a number of generated

matrices.

Keywords: parallel non-convex optimization, matrix factorization, sparse factor analysis

1 Introduction

Educational systems have witnessed a substantial transition from traditional educational methods

mainly using text books, lectures, etc. to newly developed systems which are artificial intelligent-
based systems and personally tailored to the learners [4]. Personalized Learning Systems (PLSs) and

Intelligent Tutoring Systems (ITSs) are two more well-known instances of such recently developed

educational systems. PLSs take into account learners’ individual characteristics then customize

the learning experience to the learners’ current situation and needs [2]. As computerized learning

environments, ITSs model and track student learning states [1, 6, 7]. Latent Factor Model and

Bayesian Knowledge Tracing are main classes in ITSs [3]. These new approaches encompass

computational models from different disciplines including cognitive and learning sciences, education,

computational linguistics, artificial intelligence, operations research, and other fields. More details

can be found in [1, 4–6].

Recently, [4] developed a new machine learning-based model for learning analytics, which approximate

a students knowledge of the concepts underlying a domain, and content analytics, which estimate

the relationships among a collection of questions and those concepts. This model calculates the

probability that a learner provides the correct response to a question in terms of three factors: their

understanding of a set of underlying concepts, the concepts involved in each question, and each

questions intrinsic difficulty [4]. They proposed a bi-convex maximum-likelihood-based solution to

the resulting SPARse Factor Analysis (SPARFA) problem. However, the scalability of SPARFA

when the number of questions and students significantly increase has not been studied yet.



Problem :

Sparse Matrix Factorization (SMF) is a key component in many machine learning problems

and there exist a verity a applications in real-world problems such as recommendation systems,

estimating missing values, gene expression modeling, intelligent tutoring systems (ITSs), etc.

There are different approaches to tackle with SMF rooted in linear algebra and probability

theory. In this project, given an incomplete binary matrix of students’ performances over a

set of questions, estimating the probability of success or fail over unanswered questions is of

interest. This problem is formulated using Maximum Likelihood Estimation (MLE) which leads

to a biconvex optimization problem (this formulation is based on SPARFA [4]). The resulting

optimization problem is a hard problem to deal with due to the existence of many local minima.

On the other hand, when the size of the matrix of students’ performances increase, the existing

algorithms are not successful; therefore, an efficient algorithm is required to solve this problem

for large matrices. In this project, a parallel algorithm (i.e., a parallel version of SPARFA) is

developed to solve the biconvex optimization problem and tested via a number of generated

matrices.

Keywords: parallel non-convex optimization, matrix factorization, sparse factor analysis

1 Introduction

Educational systems have witnessed a substantial transition from traditional educational methods

mainly using text books, lectures, etc. to newly developed systems which are artificial intelligent-
based systems and personally tailored to the learners [4]. Personalized Learning Systems (PLSs) and

Intelligent Tutoring Systems (ITSs) are two more well-known instances of such recently developed

educational systems. PLSs take into account learners’ individual characteristics then customize

the learning experience to the learners’ current situation and needs [2]. As computerized learning

environments, ITSs model and track student learning states [1, 6, 7]. Latent Factor Model and

Bayesian Knowledge Tracing are main classes in ITSs [3]. These new approaches encompass

computational models from different disciplines including cognitive and learning sciences, education,

computational linguistics, artificial intelligence, operations research, and other fields. More details

can be found in [1, 4–6].

Recently, [4] developed a new machine learning-based model for learning analytics, which approximate

a students knowledge of the concepts underlying a domain, and content analytics, which estimate

the relationships among a collection of questions and those concepts. This model calculates the

probability that a learner provides the correct response to a question in terms of three factors: their

understanding of a set of underlying concepts, the concepts involved in each question, and each

questions intrinsic difficulty [4]. They proposed a bi-convex maximum-likelihood-based solution to

the resulting SPARse Factor Analysis (SPARFA) problem. However, the scalability of SPARFA

when the number of questions and students significantly increase has not been studied yet.



2 Problem Statement

Let Y denote binary-valued data set of performance of N students on Q questions; therefore Y is a

matrix of size Q × N with entry Yij = 1(0) if student j answers question i correctly. Matrix Y is

highly sparse because there are a lot of unanswered questions resulted in incomplete data. One way

to estimate missing values in Y is to factorize Y into matrices W, C and M such that a function of

WC+M can estimate values of Y. It is assumed that the collection of questions is related to a small

number of abstract concepts represented by W where weight Wik (∀i = 1, . . . , Q and k = 1, . . . , K)

indicates the degree to which question i involves concept k and K is the number of latent abstract

concept. Let Ckj (∀k = 1, . . . , K and j = 1, . . . , N) denote student j’s knowledge of concept k (C is

the matrix version of Ckj ). M is an Q × N matrix representing intrinsic difficulty of each question.

It is assumed that K  Q, N so W becomes a tall, narrow Q × K matrix and C will be a short,

wide K × N matrix.






References

[1] Arthur C Graesser, Mark W Conley, and Andrew Olney, Intelligent tutoring systems., (2012).

[2] Sabine Graf et al., Personalized learning systems, Encyclopedia of the Sciences of Learning,

Springer, 2012, pp. 2594–2596.

[3] M Khajah, Rowan M Wing, Robert V Lindsey, and Michael C Mozer, Integrating latent-factor

and knowledge-tracing models to predict individual differences in learning, Proceedings of the

Seventh International Conference on Educational Data Mining, 2014.

3

[4] Andrew S Lan, Andrew E Waters, Christoph Studer, and Richard G Baraniuk, Sparse factor

analysis for learning and content analytics, The Journal of Machine Learning Research 15 (2014),

no. 1, 1959–2008.

[5] Nan Li, William W Cohen, Kenneth R Koedinger, Noboru Matsuda, and Carnegie Mellon, A

machine learning approach for automatic student model discovery., EDM, 2011, pp. 31–40.

[6] Anna N Rafferty, Emma Brunskill, Thomas L Griffiths, and Patrick Shafto, Faster teaching by

pomdp planning, Artificial intelligence in education, Springer, 2011, pp. 280–287.

[7] Shaghayegh Sahebi, Yun Huang, and Peter Brusilovsky, Predicting student performance in

solving parameterized exercises, Intelligent Tutoring Systems, Springer, 2014, pp. 496–503.
