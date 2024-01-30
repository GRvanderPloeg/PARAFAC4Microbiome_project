function result=generateNormalRandomValues(number, sigma, mu)
result = sigma.*randn(number,1) + mu;