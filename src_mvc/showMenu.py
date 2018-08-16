# This class only contains the main menu which is shown to the user.
# It has only one function the menu() function.

class showMenu:

	def menu():
		print('Choose one the following numbers:')
		print('---------------------------------')
		print('1. Read GenBank files from folder.')
		print('2. Load GenBank files.')
		print('---------------------------------')
		print('3. Run NMF model.')
		print('4. Run LDA original model.')
		print('5. Run LDA binary model.')