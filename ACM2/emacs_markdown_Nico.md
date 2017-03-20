## Basic Navigation

Begin by opening emacs. Do this by opening a terminal on the virtual machine by pressing ctrl-alt-t, then type 'emacs &' and press enter.

Once in emacs open the file named 'emacs_tutorial.txt'. The emacs mini-buffer acts like a terminal therefore you need to enter the file path and file name. You can also use tab complete.

The emacs_tutorial.txt file contains a list of commands for basic navigation of the buffers. It also contains instructions to create new buffers, move between buffers and delete buffers.

## Theme Changing

The theme of emacs can be changed to suit the users preferences. This is done by using alt-x load-theme. The options can be viewed by pressing tab. Select favoured theme then press enter.

## Calendar

A useful function of emacs is the built in calendar. To access the calendar type alt-x calendar. This will display the calendar in a new buffer. The calendar is read only however events can be added. To do this, move to the date you wish to create an event. Press i d to add event to that day, i w creates a weekly event and i y creates a yearly event. Ensure the buffer is saved. The event is saved to the buffer named 'diary', which is stored within the emacs directories.

Once the event has been created the user can then view the events of the day by typing 'alt-x diary' whilst working on another buffer.

## Tasks
A series of short tasks has been created to get the user to use emacs and learn the initial commands

1. Create new file
2. Import the ‘sample.txt’ file
3. Navigate through it
4. Add something at the end
5. Save
6. Add event to todays date
7. Go back to the file
8. Check what events are on today
9. Close the event buffer

### Solutions

1. C-x C-f then type the name of the file you want to create
2. C-x i then find the file named 'sample.txt'
3. Get used to the C-f C-b etc commands from the tutorial
4. Make an addition to the buffer
5. C-x C-s
6. M-x calendar, go to todays date, i d to insert an event to that date. Ensure the buffer is saved before exiting
7. C-x C-f find the file.
8. M-x diary whilst in the created file buffer.
9. C-o to go to the event buffer, C-x 0 to close the buffer.
