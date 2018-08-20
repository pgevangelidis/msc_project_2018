# This class only contains the main menu which is shown to the user.
# It has only one function the menu() function.

class showMenu:

	def __init__(self):
		self.choices = 8

	def menu(self):
		print('Choose one the following numbers:')
		print('----------Save/Load--------------')
		print('1. Read GenBank files from folder.')
		print('2. Load GenBank files.')
		print('---------------------------------')
		print('------------MODELS---------------')
		print('3. Run NMF model.')
		print('4. Run LDA original model.')
		print('5. Run LDA binary model.')
		print('6. Run Mixture Model')
		print('---------------------------------')
		print('-------------Plot----------------')
		print('7. Plot the lower bound curve.')
		print('---------------------------------')
		print('-------------Plot----------------')
		print('8. Exit')