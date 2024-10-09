# Setting up VS Code for Hydra Connection

First step is to install VS Code, if you do not already have it. You can go to https://code.visualstudio.com/ and click on the big Download for [Windows/Mac/Linux] button.

Once it is installed. Click on the Extension button on the left side.

![VS Code Extension Button](screenshots/extension_button.png)

Then search for the word "remote" to pull up a bunch of results. Click Install on the first option, which should be "Remote Development", which is an extension pack that will install multiple extensions.

![VS Code Extension Search Remote](screenshots/remote_search.png)

If this worked, it will show that multiple extensions are now installed.

![VS Code Remote Extensions Installed](screenshots/remote_installed.png)

Now click on the new Remote Explorer icon on the left sidebar. Then click on the + icon next to SSH to add an SSH connection to Hydra.

![VS Code Remote Panel Blank](screenshots/remote_panel.png)

A little panel will show up in the top middle of VS Code with an example SSH connection command. Fill it in with `ssh USER@hydra-login01.si.edu`.

![VS Code Add SSH](screenshots/add_ssh.png)

You will then get a follow-up prompt that asks where to store this config. Just choose whatever is the first option.

![VS Code Add SSH Config](screenshots/add_ssh_config.png)

It will do a little bit of thinking, but then you should see the Hydra SSH connection added to your Remote Explorer tab.

Click on the icon either for "Connect in Current Window" or "Connect in New Window" to connect to Hydra. It doesn't really matter which one you choose.

![VS Code Remote Panel with Hydra](screenshots/remote_panel_hydra.png)

It will again open up a little prompt in the top middle of VS Code, which asks for your Hydra password. Enter it here, and press enter.

![VS Code Hydra Password Prompt](screenshots/hydra_password_prompt.png)

If this is the first time using VS Code on Hydra, it will take a while to install VS Code helper files on Hydra in your space. But when it finishes, you should see the options to Open Folder, and it will also open up a Terminal window for you.

![VS Code Remote Open Folder](screenshots/remote_open_folder.png)

Here is an example screenshot, after opening my own Home folder. I have clicked on a random file that shows how it opens it up in the Editor pane. I also have entered a command in the terminal, that shows that the path between the Explorer pane and the terminal location can get slightly disconnected.

![VS Code Remote Full Window](screenshots/remote_full_window.png)
