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

1

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

1

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
