# NeuralNet subfcn

#define list of crossvalidation flags that are callable
FLAGS = flags(
	flag_numeric('width1', 16),
	flag_numeric('dropout', .1),
	flag_numeric('depth', 1),
	flag_numeric('epoch_num', 64),
	flag_numeric('batch_size', 10),
	flag_string('activationfcn', 'relu')
)

#sequentially build model

if(FLAGS$depth <= 1){
	model = keras_model_sequential() %>% layer_dense(units = 1, activation = 'linear', input_shape = c(num_zs))
} else{
	model = keras_model_sequential() %>% layer_dense(units = FLAGS$width1, input_shape = c(num_zs)) %>%
	layer_activation_leaky_relu() %>%
	layer_dropout(rate = FLAGS$dropout)
	if(FLAGS$depth > 2){
for(laycount in 1:(FLAGS$depth - 2)){
	model %>% layer_dense(units = FLAGS$width1) %>%
	layer_activation_leaky_relu() %>%
	layer_dropout(rate = FLAGS$dropout)
	}
}
model %>% layer_dense(units = 1, activation = 'linear')}



compile(model, optimizer = 'adam', loss = 'mse')

history = model %>% fit(x = train_dat, y =train_lab, epochs = FLAGS$epoch_num, batch_size = FLAGS$batch_size, verbose = 0)

evaluate(model, test_dat, test_lab, verbose = 0)

k_clear_session()