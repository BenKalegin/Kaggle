install.packages("drat", repos = "https://cran.rstudio.com")
drat:::addRepo("dmlc")
install.packages("mxnet")
require(mxnet)


data <- mx.symbol.Variable('data')
# 1st convolutional layer
conv_1 <- mx.symbol.Convolution(data = data, kernel = c(3), num_filter = 16)
pool_1 <- mx.symbol.Pooling(data = conv_1, pool_type = "max", kernel = c(2))
flatten <- mx.symbol.Flatten(data = pool_1)
dotplus1 <- mx.symbol. dot(flatten, )



tanh_1 <- mx.symbol.Activation(data = conv_1, act_type = "tanh")
pool_1 <- mx.symbol.Pooling(data = tanh_1, pool_type = "max", kernel = c(2, 2), stride = c(2, 2))
# 2nd convolutional layer
conv_2 <- mx.symbol.Convolution(data = pool_1, kernel = c(5, 5), num_filter = 50)
tanh_2 <- mx.symbol.Activation(data = conv_2, act_type = "tanh")
pool_2 <- mx.symbol.Pooling(data = tanh_2, pool_type = "max", kernel = c(2, 2), stride = c(2, 2))
# 1st fully connected layer
flatten <- mx.symbol.Flatten(data = pool_2)
fc_1 <- mx.symbol.FullyConnected(data = flatten, num_hidden = 500)
tanh_3 <- mx.symbol.Activation(data = fc_1, act_type = "tanh")
# 2nd fully connected layer
fc_2 <- mx.symbol.FullyConnected(data = tanh_3, num_hidden = 40)
# Output. Softmax output since we'd like to get some probabilities.
NN_model <- mx.symbol.SoftmaxOutput(data = fc_2)

# Pre-training set up
#-------------------------------------------------------------------------------

# Set seed for reproducibility
mx.set.seed(100)

# Device used. CPU in my case.
devices <- mx.cpu()

# Training
#-------------------------------------------------------------------------------

# Train the model
model <- mx.model.FeedForward.create(NN_model,
                                     X = train_array,
                                     y = train_y,
                                     ctx = devices,
                                     num.round = 480,
                                     array.batch.size = 40,
                                     learning.rate = 0.01,
                                     momentum = 0.9,
                                     eval.metric = mx.metric.accuracy,
                                     epoch.end.callback = mx.callback.log.train.metric(100))

# Testing
#-------------------------------------------------------------------------------

# Predict labels
predicted <- predict(model, test_array)
# Assign labels
predicted_labels <- max.col(t(predicted)) - 1
# Get accuracy
sum(diag(table(test[, 1], predicted_labels))) / 40


width = 100;

height = 1;

net = NetChain[{
    ConvolutionLayer[16( * number of filters * ), { height, 3( * kernel size * ) }],
        PoolingLayer[{ height, 2( * pooling size * ) }],
            FlattenLayer[],
            DotPlusLayer[10],
            DotPlusLayer[2( * number of classes * )],
        SoftmaxLayer["Output" -> NetDecoder[{ "Class", { "tech", "banks" }}]]
    },
  "Input" -> NetEncoder[{ "Image", { width, height }, ColorSpace -> "Grayscale" }]
  ]